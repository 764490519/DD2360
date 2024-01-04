#### Results of Test：

##### **only use quick select in kernel, use data of 1024k_12, find 10 min element**

for algorithm of quick select

| Test | Copy time host to device | Kernel time 1-calculate the distance | Kernel time 2-calculate the min-k | Copy time device to Host |
| ---- | ------------------------ | ------------------------------------ | --------------------------------- | ------------------------ |
| 1    | 0.002031                 | 0.000371                             | 0.025717                          | 0.000043                 |
| 2    | 0.002355                 | 0.000352                             | 0.025666                          | 0.000039                 |
| 3    | 0.002490                 | 0.000327                             | 0.025752                          | 0.000040                 |
| 4    | 0.002455                 | 0.000314                             | 0.025758                          | 0.000041                 |

for origin 

| Test | Copy H to D | Kernel time | Copy D to H | calculate min-k |
| ---- | ----------- | ----------- | ----------- | --------------- |
| 1    | 0.002381    | 0.000477    | 0.003549    | 0.031021        |
| 2    | 0.002331    | 0.000479    | 0.003291    | 0.028545        |
| 3    | 0.002384    | 0.000487    | 0.003213    | 0.034279        |
| 4    | 0.02712     | 0.000440    | 0.004715    | 0.027563        |

##### **for this data but find 1000 min element**

for quickselect

| Test | Copy time host to device | Kernel time 1-calculate the distance | Kernel time 2-calculate the min-k | Copy time device to Host |
| ---- | ------------------------ | ------------------------------------ | --------------------------------- | ------------------------ |
| 1    | 0.002425                 | 0.000348                             | 0.158351                          | 0.000062                 |
| 2    | 0.002709                 | 0.000373                             | 0.161793                          | 0.000062                 |

for origin

| Test | Copy H to D | Kernel time | Copy D to H | calculate min-k |
| ---- | ----------- | ----------- | ----------- | --------------- |
| 1    | 0.002022    | 0.000465    | 0.003400    | 3.207194        |
| 2    | 0.002154    | 0.000454    | 0.004254    | 2.904574        |

##### for a much larger data 5000k_512 and find min 10000

Obviously, it is much faster than the code rodinia support.

for quickselect

| Test | Copy time host to device | Kernel time 1-calculate the distance | Kernel time 2-calculate the min-k | Copy time device to Host |
| ---- | ------------------------ | ------------------------------------ | --------------------------------- | ------------------------ |
| 1    | 0.019673                 | 0.000627                             | 0.951120                          | 0.000082                 |
| 2    | 0.010281                 | 0.000530                             | 0.969164                          | 0.000081                 |

for origin code

| Test | Copy H to D | Kernel time | Copy D to H | calculate min-k |
| ---- | ----------- | ----------- | ----------- | --------------- |
| 1    | 0.009263    | 0.000716    | 0.004688    | 146.812323      |
| 2    | 0.010410    | 0.000661    | 0.004430    | 143.467046      |

The difference is very obvious. On average, the quick selection algorithm is 151 times faster than the original algorithm.

#### Conclusion

The improvements I made mainly focused on the quick selection algorithm.

Firstly, make the calculation of the nearest neighbors to the kernel can reduce the time that copy distance to host. Only need to copy the result to Host.

Secondly, I consider which algorithm is appropriate to find the k-th min element. There are several options:

1. quickselect: On average it has O(n) time complexity. In kernel, I can break it into m blocks, and find k-th min element in each block. Then use a thread to calculate the k-th min element in total. Then the time complexity is : (n numbers, k-th min, m blocks,Ignore constants)

$$
O(\frac n m+mk)
$$

so it is easy to get the value of m make the time complexity minimum. 
$$
m  = \sqrt \frac n k\\
\frac n m+mk = 2\sqrt{nk}
\\So, the \ time \ complexity \ is \  O(\sqrt{nk})
$$

2. Min-Heap: Obviously, O(n log k). If still make it into m blocks, the merging heap cost mk log k, unless use Fibonacci Tree, which merge only use O(1). ease to get the minimum time complexity:

   
   $$
   O(\frac n m logk+mklog k)\\
   when \ m=\sqrt \frac n k\\
   get\ O(\sqrt{nk}logk)
   $$
   

   There are quite some difficulties to achieve heap in kernel and it is not better than method 1.

3. sort algorithm: worse than heap because the complexity of O(nlogn)

Above all, use quickselect and break distances array into m blocks.