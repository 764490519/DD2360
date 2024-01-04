 /*
 * nn.cu
 * Nearest Neighbor
 *
 */

#include <stdio.h>
#include <sys/time.h>
#include <float.h>
#include <vector>
#include "cuda.h"

#define min( a, b )			a > b ? b : a
#define ceilDiv( a, b )		( a + b - 1 ) / b
#define print( x )			printf( #x ": %lu\n", (unsigned long) x )
#define DEBUG				false

#define DEFAULT_THREADS_PER_BLOCK 256

#define MAX_ARGS 10
#define REC_LENGTH 53 // size of a record in db
#define LATITUDE_POS 28	// character position of the latitude value in each record
#define OPEN 10000	// initial value of nearest neighbors
typedef struct latLong
{
  float lat;
  float lng;
} LatLong;

typedef struct record
{
  char recString[REC_LENGTH];
  float distance;
} Record;

int loadData(char *filename,std::vector<Record> &records,std::vector<LatLong> &locations);
void findLowest(std::vector<Record> &records,float *distances,int numRecords,int topN);
void printUsage();
int parseCommandline(int argc, char *argv[], char* filename,int *r,float *lat,float *lng,
                     int *q, int *t, int *p, int *d);
void printLowest(std::vector<Record> &records, int *min_record, int topN, float *dis_min);
/**
* Kernel
* Executed on GPU
* Calculates the Euclidean distance from each record in the database to the target position
*/
__global__ void euclid(LatLong *d_locations, float *d_distances, int numRecords,float lat, float lng)
{
	//int globalId = gridDim.x * blockDim.x * blockIdx.y + blockDim.x * blockIdx.x + threadIdx.x;
	int globalId = blockDim.x * ( gridDim.x * blockIdx.y + blockIdx.x ) + threadIdx.x; // more efficient
    LatLong *latLong = d_locations+globalId;
    if (globalId < numRecords) {
        float *dist=d_distances+globalId;
        *dist = (float)sqrt((lat-latLong->lat)*(lat-latLong->lat)+(lng-latLong->lng)*(lng-latLong->lng));
	}
}

//added by yitong
__device__ void swap(float* a, float* b) {
    float t = *a;
    *a = *b;
    *b = t;
}

/*
partition the array, and return the pivot index
*/
__device__ int partition(float* arr, int low, int high) {
    float pivot = arr[high];
    int i = low;
    for (int j = low; j < high; j++) {
        if (arr[j] < pivot) {
            swap(&arr[i], &arr[j]);
            i++;
        }
    }
    swap(&arr[i], &arr[high]);
    return i;
}
/*
find the min "numMin" elements in the array from "offset_mul" to "offset_mul + offset"
the result is stored in "d_minLoc" from "offset_min" to "offset_min + numMin"
time complexity is O(numRecords/sqrt(numRecords/numMin))
find_min and find_min_final uses the same algorithm, Quickselect
*/

__global__ void find_min( float *d_distances, int numRecords, int *d_minLoc, int offset, int numMin, float *tmp)
{
    int globalId = blockIdx.x * blockDim.x + threadIdx.x;
    int low = 0, high = (globalId + 1)* offset > numRecords ? numRecords - globalId* offset - 1 : offset - 1;
    int max = high, offset_min = numMin * globalId, k = 0, offset_mul = offset * globalId, pivotIndex;
    float min_k;
    if(offset_mul < numRecords){
        float* now = tmp + offset_mul;
        for(int i = 0; i < max; i++){
            now[i] = d_distances[i + offset_mul];
        }
        //quickselect
        while (low <= high) {
            // Partition the array and get the pivot index
            pivotIndex = partition(now, low, high);
           if (pivotIndex == numMin){
                min_k = now[pivotIndex];
                break;
            }   
            // If numMin is less, continue to the left part
            else if (pivotIndex > numMin) high = pivotIndex - 1;

            // If numMin is more, continue to the right part
            else low = pivotIndex + 1;
        }
        //get the min "numMin" elements
        for(int i = 0; i < max; i++){
            if(d_distances[i+offset_mul] < min_k){
                d_minLoc[offset_min + k] = i + offset_mul;
                k++;
            }
        }
        //if there are elements equal in min_k, get them 
        if(k < numMin){
            for(int i = 0; i < max; i++){
                if((d_distances[i+offset_mul] == min_k) && k < numMin){
                    d_minLoc[offset_min + k] = i + offset_mul;
                    k++;
                }
            }
        }      
    }
  //printf("\n");
}
/*
find the min "numMin" element on the base of the result of "find_min"
the result is stored in "d_minLoc" from "0" to "numMin"
time complexity is O(sqrt(numRecords/numMin) * numMin)
*/
__global__ void find_min_final(float *d_distances, int num, int *d_minLoc, int numMin, float *min_dis, float *dis_min, int *d_minmem)
{
    float min_k = 0;
    int pivotIndex;
    int low = 0, high = num - 1, k = 0;
    for(int i = 0; i < num; i++)
    {
        dis_min[i] = d_distances[d_minLoc[i]];
    }
    //quickselect
    while (low <= high) {
        // Partition the array and get the pivot index
        pivotIndex = partition(dis_min, low, high);
         // If pivot itself is the kth smallest element
        if (pivotIndex == numMin){
            min_k = dis_min[pivotIndex];
            break;
        }   
         // If k is less, continue to the left part
        else if (pivotIndex > numMin) high = pivotIndex - 1;

         // If k is more, continue to the right part
        else low = pivotIndex + 1;
    }
    //get the min "numMin" elements
    for(int i = 0; i < num; i++)
    {
        if(d_distances[d_minLoc[i]] < min_k)
        {
          d_minmem[k] = d_minLoc[i];
          printf("d_minmem[%d] :%d\n", i,d_minmem[k]);
          k++;
        }
    }
    //if there are elements equal in min_k, get them
    if(k < numMin){
        for(int i = 0; i < num; i++){
          if((d_distances[d_minLoc[i]] == min_k )&& (k < numMin)){
            d_minmem[k] = d_minLoc[i];
          printf("d_minmem[%d] :%d\n", i,d_minmem[k]);
            k++;
          }
        }
      }     
    for(int i = 0; i< numMin; i++){
      min_dis[i] = d_distances[d_minmem[i]];
      printf("min_dis[%d] :%f\n", i,min_dis[i]);
    }
}
double cpuSecond() {
   cudaDeviceSynchronize();
   struct timeval tp;
   gettimeofday(&tp,NULL);
   return ((double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}

/**
* This program finds the k-nearest neighbors
**/

int main(int argc, char* argv[])
{
	float lat, lng;
	int quiet=0,timing=0,platform=0,device=0;

  std::vector<Record> records;
	std::vector<LatLong> locations;
	char filename[100];
	int resultsCount=10;
  float *dis_min , *tmp;

    // parse command line
    if (parseCommandline(argc, argv, filename,&resultsCount,&lat,&lng,
                     &quiet, &timing, &platform, &device)) {
      printUsage();
      return 0;
    }

    int numRecords = loadData(filename,records,locations);
    if (resultsCount > numRecords) resultsCount = numRecords;

    //for(i=0;i<numRecords;i++)
    //  printf("%s, %f, %f\n",(records[i].recString),locations[i].lat,locations[i].lng);


    //Pointers to host memory
	//float *distances;
	//Pointers to device memory
	LatLong *d_locations;
	float *d_distances, *min_record_dis, *tmp_final_float;
  int *d_minLoc, *min_record, *d_minmem;


	// Scaling calculations - added by Sam Kauffman
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties( &deviceProp, 0 );
	cudaThreadSynchronize();
	unsigned long maxGridX = deviceProp.maxGridSize[0];
	unsigned long threadsPerBlock = min( deviceProp.maxThreadsPerBlock, DEFAULT_THREADS_PER_BLOCK );
	size_t totalDeviceMemory;
	size_t freeDeviceMemory;
	cudaMemGetInfo(  &freeDeviceMemory, &totalDeviceMemory );
	cudaThreadSynchronize();
	unsigned long usableDeviceMemory = freeDeviceMemory * 85 / 100; // 85% arbitrary throttle to compensate for known CUDA bug
	unsigned long maxThreads = usableDeviceMemory / 12; // 4 bytes in 3 vectors per thread
	if ( numRecords > maxThreads )
	{
		fprintf( stderr, "Error: Input too large.\n" );
		exit( 1 );
	}
	unsigned long blocks = ceilDiv( numRecords, threadsPerBlock ); // extra threads will do nothing
	unsigned long gridY = ceilDiv( blocks, maxGridX );
	unsigned long gridX = ceilDiv( blocks, gridY );
	// There will be no more than (gridY - 1) extra blocks
	dim3 gridDim( gridX, gridY );
  dim3 grid_min(((int)(sqrt((float)numRecords/(float)resultsCount)+1) + threadsPerBlock - 1)/threadsPerBlock);
	if ( DEBUG )
	{
		print( totalDeviceMemory ); // 804454400
		print( freeDeviceMemory );
		print( usableDeviceMemory );
		print( maxGridX ); // 65535
		print( deviceProp.maxThreadsPerBlock ); // 1024
		print( threadsPerBlock );
		print( maxThreads );
		print( blocks ); // 130933
		print( gridY );
		print( gridX );
	}

	/**
	* Allocate memory on host and device
	*/
	//distances = (float *)malloc(sizeof(float) * numRecords);
  min_record = (int *)malloc(sizeof(int) * resultsCount);
  dis_min = (float *)malloc(sizeof(float) * resultsCount);
	cudaMalloc((void **) &d_locations,sizeof(LatLong) * numRecords);
	cudaMalloc((void **) &d_distances,sizeof(float) * numRecords);
  cudaMalloc((void **) &d_minLoc,sizeof(int) * (int)(sqrt((float)numRecords/(float)resultsCount)+1)* resultsCount);
  cudaMalloc((void **) &min_record_dis,sizeof(float) * resultsCount);
  cudaMalloc((void **) &tmp,sizeof(float) * numRecords);
  cudaMalloc((void **) &tmp_final_float,sizeof(float) * (int)(sqrt((float)numRecords/(float)resultsCount)+1)* resultsCount);
  cudaMalloc((void **) &d_minmem,sizeof(int) * resultsCount);
   /**
    * Transfer data from host to device
    */
    cudaMemcpy( d_locations, &locations[0], sizeof(LatLong) * numRecords, cudaMemcpyHostToDevice);
    /**
    * Execute kernel
    */
    euclid<<< gridDim, threadsPerBlock >>>(d_locations,d_distances,numRecords,lat,lng);
    cudaThreadSynchronize();
    int offset =(float)numRecords / (int)(sqrt((float)numRecords/(float)resultsCount)+1) + 1;
    find_min<<< grid_min , threadsPerBlock >>>(d_distances, numRecords, d_minLoc, offset, resultsCount,tmp);
    cudaThreadSynchronize();
    find_min_final<<< 1, 1 >>>(d_distances, (int)(sqrt((float)numRecords/(float)resultsCount)+1) * resultsCount, d_minLoc, resultsCount, min_record_dis, tmp_final_float, d_minmem);
    cudaThreadSynchronize();
    //Copy data from device memory to host memory
    //cudaMemcpy( distances, d_distances, sizeof(float)*numRecords, cudaMemcpyDeviceToHost );
    cudaMemcpy( min_record, d_minmem, sizeof(int)*resultsCount, cudaMemcpyDeviceToHost );
    cudaMemcpy( dis_min, min_record_dis, sizeof(float)*resultsCount, cudaMemcpyDeviceToHost );
	// find the resultsCount least distances
    //findLowest(records,distances,numRecords,resultsCount);
    printLowest(records, min_record, resultsCount, dis_min);
    // print out results
    /*if (!quiet)
    for(i=0;i<resultsCount;i++) {
      printf("%s --> Distance=%f\n",records[i].recString,records[i].distance);
    }*/
    free(min_record);
    free(dis_min);
    //Free memory
	cudaFree(d_locations);
	cudaFree(d_distances);
  cudaFree(d_minLoc);
  cudaFree(min_record_dis);
  cudaFree(tmp);
  cudaFree(tmp_final_float);
  cudaFree(d_minmem);
    return 0;

}

int loadData(char *filename,std::vector<Record> &records,std::vector<LatLong> &locations){
    FILE   *flist,*fp;
	int    i=0;
	char dbname[64];
	int recNum=0;

    /**Main processing **/

    flist = fopen(filename, "r");
	while(!feof(flist)) {
		/**
		* Read in all records of length REC_LENGTH
		* If this is the last file in the filelist, then done
		* else open next file to be read next iteration
		*/
		if(fscanf(flist, "%s\n", dbname) != 1) {
            fprintf(stderr, "error reading filelist\n");
            exit(0);
        }
        fp = fopen(dbname, "r");
        if(!fp) {
            printf("error opening a db\n");
            exit(1);
        }
        // read each record
        while(!feof(fp)){
            Record record;
            LatLong latLong;
            fgets(record.recString,49,fp);
            fgetc(fp); // newline
            if (feof(fp)) break;

            // parse for lat and long
            char substr[6];

            for(i=0;i<5;i++) substr[i] = *(record.recString+i+28);
            substr[5] = '\0';
            latLong.lat = atof(substr);

            for(i=0;i<5;i++) substr[i] = *(record.recString+i+33);
            substr[5] = '\0';
            latLong.lng = atof(substr);

            locations.push_back(latLong);
            records.push_back(record);
            recNum++;
        }
        fclose(fp);
    }
    fclose(flist);
//    for(i=0;i<rec_count*REC_LENGTH;i++) printf("%c",sandbox[i]);
    return recNum;
}

void findLowest(std::vector<Record> &records,float *distances,int numRecords,int topN){
  int i,j;
  float val;
  int minLoc;
  Record *tempRec;
  float tempDist;

  for(i=0;i<topN;i++) {
    minLoc = i;
    for(j=i;j<numRecords;j++) {
      val = distances[j];
      if (val < distances[minLoc]) minLoc = j;
    }
    // swap locations and distances
    tempRec = &records[i];
    records[i] = records[minLoc];
    records[minLoc] = *tempRec;

    tempDist = distances[i];
    distances[i] = distances[minLoc];
    distances[minLoc] = tempDist;

    // add distance to the min we just found
    records[i].distance = distances[i];
  }
}
void printLowest(std::vector<Record> &records, int *min_record, int topN, float *dis_min){
  for(int i = 0; i < topN; i++)
    printf("%s --> Distance=%f\n",records[min_record[i]].recString,dis_min[i]);
}
int parseCommandline(int argc, char *argv[], char* filename,int *r,float *lat,float *lng,
                     int *q, int *t, int *p, int *d){
    int i;
    if (argc < 2) return 1; // error
    strncpy(filename,argv[1],100);
    char flag;

    for(i=1;i<argc;i++) {
      if (argv[i][0]=='-') {// flag
        flag = argv[i][1];
          switch (flag) { 
            case 'r': // number of results
              i++;
              *r = atoi(argv[i]);
              break;
            case 'l': // lat or lng
              if (argv[i][2]=='a') {//lat
                *lat = atof(argv[i+1]);
              }
              else {//lng
                *lng = atof(argv[i+1]);
              }
              i++;
              break;
            case 'h': // help
              return 1;
            case 'q': // quiet
              *q = 1;
              break;
            case 't': // timing
              *t = 1;
              break;
            case 'p': // platform
              i++;
              *p = atoi(argv[i]);
              break;
            case 'd': // device
              i++;
              *d = atoi(argv[i]);
              break;
        }
      }
    }
    if ((*d >= 0 && *p<0) || (*p>=0 && *d<0)) // both p and d must be specified if either are specified
      return 1;
    return 0;
}

void printUsage(){
  printf("Nearest Neighbor Usage\n");
  printf("\n");
  printf("nearestNeighbor [filename] -r [int] -lat [float] -lng [float] [-hqt] [-p [int] -d [int]]\n");
  printf("\n");
  printf("example:\n");
  printf("$ ./nearestNeighbor filelist.txt -r 5 -lat 30 -lng 90\n");
  printf("\n");
  printf("filename     the filename that lists the data input files\n");
  printf("-r [int]     the number of records to return (default: 10)\n");
  printf("-lat [float] the latitude for nearest neighbors (default: 0)\n");
  printf("-lng [float] the longitude for nearest neighbors (default: 0)\n");
  printf("\n");
  printf("-h, --help   Display the help file\n");
  printf("-q           Quiet mode. Suppress all text output.\n");
  printf("-t           Print timing information.\n");
  printf("\n");
  printf("-p [int]     Choose the platform (must choose both platform and device)\n");
  printf("-d [int]     Choose the device (must choose both platform and device)\n");
  printf("\n");
  printf("\n");
  printf("Notes: 1. The filename is required as the first parameter.\n");
  printf("       2. If you declare either the device or the platform,\n");
  printf("          you must declare both.\n\n");
}
