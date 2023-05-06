package math;


public final class Sort
{

 public Sort()
 {
 }

 static int pivot(double a[], int i, int j)
 {
     int k;
     for(k = i + 1; k <= j && a[i] == a[k]; k++);
     if(k > j)
         return -1;
     if(a[i] >= a[k])
         return i;
     else
         return k;
 }

 static int pivot(int a[], int i, int j)
 {
     int k;
     for(k = i + 1; k <= j && a[i] == a[k]; k++);
     if(k > j)
         return -1;
     if(a[i] >= a[k])
         return i;
     else
         return k;
 }

 static int partition(double a[], int sort[], int i, int j, double x)
 {
     int l = i;
     int r = j;
     do
     {
         if(l > r)
             break;
         for(; l <= j && a[l] < x; l++);
         for(; r >= i && a[r] >= x; r--);
         if(l > r)
             break;
         int ts = sort[l];
         double t = a[l];
         sort[l] = sort[r];
         a[l] = a[r];
         sort[r] = ts;
         a[r] = t;
         l++;
         r--;
     } while(true);
     return l;
 }

 static int partition(double a[], int i, int j, double x)
 {
     int l = i;
     int r = j;
     do
     {
         if(l > r)
             break;
         for(; l <= j && a[l] < x; l++);
         for(; r >= i && a[r] >= x; r--);
         if(l > r)
             break;
         double t = a[l];
         a[l] = a[r];
         a[r] = t;
         l++;
         r--;
     } while(true);
     return l;
 }

 static int partition(int a[], int sort[], int i, int j, int x)
 {
     int l = i;
     int r = j;
     do
     {
         if(l > r)
             break;
         for(; l <= j && a[l] < x; l++);
         for(; r >= i && a[r] >= x; r--);
         if(l > r)
             break;
         int ts = sort[l];
         int t = a[l];
         sort[l] = sort[r];
         a[l] = a[r];
         sort[r] = ts;
         a[r] = t;
         l++;
         r--;
     } while(true);
     return l;
 }

 static int partition(int a[], int i, int j, int x)
 {
     int l = i;
     int r = j;
     do
     {
         if(l > r)
             break;
         for(; l <= j && a[l] < x; l++);
         for(; r >= i && a[r] >= x; r--);
         if(l > r)
             break;
         int t = a[l];
         a[l] = a[r];
         a[r] = t;
         l++;
         r--;
     } while(true);
     return l;
 }

 public static void quick(double a[], int sort[], int i, int j)
 {
     if(i == j)
         return;
     int p = pivot(a, i, j);
     if(p != -1)
     {
         int k = partition(a, sort, i, j, a[p]);
         quick(a, sort, i, k - 1);
         quick(a, sort, k, j);
     }
 }

 public static void quick(double a[], int i, int j)
 {
     if(i == j)
         return;
     int p = pivot(a, i, j);
     if(p != -1)
     {
         int k = partition(a, i, j, a[p]);
         quick(a, i, k - 1);
         quick(a, k, j);
     }
 }

 public static void quick(int a[], int sort[], int i, int j)
 {
     if(i == j)
         return;
     int p = pivot(a, i, j);
     if(p != -1)
     {
         int k = partition(a, sort, i, j, a[p]);
         quick(a, sort, i, k - 1);
         quick(a, sort, k, j);
     }
 }

 public static void quick(int a[], int i, int j)
 {
     if(i == j)
         return;
     int p = pivot(a, i, j);
     if(p != -1)
     {
         int k = partition(a, i, j, a[p]);
         quick(a, i, k - 1);
         quick(a, k, j);
     }
 }

 public static void quick(double a[], int sort[])
 {
     quick(a, sort, 0, a.length - 1);
 }

 public static void quick(double a[])
 {
     quick(a, 0, a.length - 1);
 }

 public static void quick(int a[], int sort[])
 {
     quick(a, sort, 0, a.length - 1);
 }

 public static void quick(int a[])
 {
     quick(a, 0, a.length - 1);
 }

 public static void bubble()
 {
 }

 public static void nanchatte()
 {
 }
}
