#define _CRT_SECURE_NO_WARNINGS 1
#include <iostream>
#include<string>
#include<fstream>
#include <time.h>
#include<stdlib.h>
#include<math.h>
#include <cstring>

using namespace std;

int a[500001] = { 0 };
int n = 0;
unsigned long long int cnt_compare = 0;
int cnt_assign = 0;
clock_t begin;
//stable sorting: bubble sort, insertion sort, merge sort
bool isNum(string s)
{
    for (int i = 0; i < s.length(); i++)
    {
        if (s[i] < '0' || s[i] > '9')
            return false;
    }
    return true;
}
void write_input(int a[], int size, int datatype)
{
    FILE* f = NULL;
    if (datatype == 1)f = fopen("input_1.txt", "w");
    if (datatype == 2)f = fopen("input_2.txt", "w");
    if (datatype == 3)f = fopen("input_3.txt", "w");
    if (datatype == 4)f = fopen("input_4.txt", "w");

    for (int i = 0; i <= size; i++)
    {
        fprintf(f, "%i", a[i]);
        fprintf(f, " ");
    }
    fclose(f);
}

template <class T>
void HoanVi(T& a, T& b)
{
    T x = a;
    a = b;
    b = x;
}

void GenerateRandomData(int a[], int n)
{
    srand((unsigned int)time(NULL));
    for (int i = 0; i < n; i++)
    {
        a[i] = rand() % n;
    }
}

void GenerateSortedData(int a[], int n)
{
    for (int i = 0; i < n; i++)
    {
        a[i] = i;
    }
}

void GenerateReverseData(int a[], int n)
{
    for (int i = 0; i < n; i++)
    {
        a[i] = n - 1 - i;
    }
}

void GenerateNearlySortedData(int a[], int n)
{
    for (int i = 0; i < n; i++)
    {
        a[i] = i;
    }
    srand((unsigned int)time(NULL));
    for (int i = 0; i < 10; i++)
    {
        int r1 = rand() % n;
        int r2 = rand() % n;
        while (r1 == r2)
        {
            int r1 = rand() % n;
            int r2 = rand() % n;
        }
        HoanVi(a[r1], a[r2]);
    }
}

void GenerateData(int a[], int n, int dataType)
{
    switch (dataType)
    {
    case 0:	// ngau nhien
        GenerateRandomData(a, n);
        break;
    case 1:	// co thu tu
        GenerateSortedData(a, n);
        break;
    case 2:	// co thu tu nguoc
        GenerateReverseData(a, n);
        break;
    case 3:	// gan nhu sort
        GenerateNearlySortedData(a, n);
        break;
    default:
        printf("Error: unknown data type!\n");
    }
}

void Heapify(int a[], int n, int i)
{
    int largest = i;
    int left = i * 2 + 1;
    int right = i * 2 + 2;
    if (++cnt_compare && left < n && ++cnt_compare && a[i] > a[left])
    {
        largest = left;
    }
    if (++cnt_compare && right < n && ++cnt_compare && a[right] < a[largest])
    {
        largest = right;
    }
    if (++cnt_compare && largest != i)
    {
        swap(a[i], a[largest]);
        Heapify(a, n, largest);
    }
}
void buildMinHeap(int a[], int n)
{
    for (int i = (n / 2) - 1; ++cnt_compare, i >= 0; i--)
    {
        Heapify(a, n, i);
    }
}
void heapSort(int a[], int n)
{
    buildMinHeap(a, n);
    /*for (int i = 0; i < n; i++)
    {
        cout << a[i] << ' ';
    }
    cout << endl << 1 << endl;*/

    for (int i = n - 1; ++cnt_compare, i > 0; i--)
    {
        swap(a[0], a[i]);
        Heapify(a, i, 0);
    }
}

void shakerSort(int a[], int n)
{
    int left = 0;
    int right = n - 1;
    int k = 0;
    while (++cnt_compare && left < right)
    {
        for (int i = left; cnt_compare++,i < right; i++)
        {
            if (++cnt_compare && a[i] > a[i + 1])
            {
                k = i;
                swap(a[i], a[i + 1]);
            }
        }
        right = k;
        for (int i = right; cnt_compare++, i > left; i--)
        {
            if (++cnt_compare && a[i] < a[i - 1])
            {
                k = i;
                swap(a[i], a[i - 1]);
            }
        }
        left = k;
    }
}

void bubbleSort(int a[], int n)
{
    for (int i = n - 1; i >= 0; i--)
    {
        bool swapped = false;
        //day max ve cuoi mang
        for (int j = 0; j < i; j++)
            if (++cnt_compare && a[j] > a[j + 1])
            {
                swapped = true;
                HoanVi(a[j], a[j + 1]);
            }
        if (!swapped)
        {
            return;
        }
    }
}
int partition(int a[], int first, int last)
{
    int pivot = a[(first + last) / 2];
    int i = first;
    int j = last;
    int tmp;
    while (++cnt_compare && i <= j)
    {
        while (++cnt_compare && a[i] < pivot)
        {
            i++;
        }
        while (++cnt_compare && a[j] > pivot)
        {
            j--;
        }
        if (++cnt_compare && i <= j)
        {
            tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
            i++;
            j--;
        }
    }
    return i;
}


void QuickSort(int a[], int low, int high)
{
    if (++cnt_compare && low < high)
    {
        int p = partition(a, low, high);
        QuickSort(a, low, p - 1);
        QuickSort(a, p + 1, high);

    }
}


void quickSort(int a[], int n)
{
    QuickSort(a, 0, n - 1);
}

void countingSort(int a[], int n)
{
    int max = a[0];
    for (int i = 1; ++cnt_compare && i < n; i++)
    {
        if (++cnt_compare && a[i] > max)
        {
            max = a[i];
        }
    }

    int* count = new int[max + 1];
    for (int i = 0; ++cnt_compare && i <= max; i++)
    {
        count[i] = 0;
    }

    for (int i = 0; ++cnt_compare && i < n; i++)
    {
        count[a[i]]++;
    }

    for (int i = 1; ++cnt_compare && i <= max; i++)
    {
        count[i] += count[i - 1];
    }

    int* temp = new int[n];
    for (int i = 0; ++cnt_compare && i < n; i++)
    {
        temp[count[a[i]] - 1] = a[i];
        count[a[i]]--;
    }

    for (int i = 0; ++cnt_compare && i < n; i++)
    {
        a[i] = temp[i];
    }
    delete[] count;
    delete[] temp;
}

void insertionSort(int arr[], int n)	////ham duoc tham khao tu geeksforgeeks.org
{
    cnt_compare = 0;
    int i, key, j;
    for (i = 1; ++cnt_compare && i < n; i++)
    {
        key = arr[i];
        j = i - 1;
        // Move elements of arr[0..i-1], 
      // that are greater than key, to one
      // position ahead of their
      // current position
        while (++cnt_compare && j >= 0 && ++cnt_compare && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

void shellSort(int arr[], int n)	//ham duoc tham khao tu geeksforgeeks.org
{
    // Start with a big gap, then reduce the gap
    for (int gap = n / 2; ++cnt_compare && gap > 0; gap /= 2)
    {
        //Do a gapped insertion sort for this gap size.
            // The first gap elements a[0..gap-1] are already in gapped order
            // keep adding one more element until the entire array is
            // gap sorted 
        for (int i = gap; ++cnt_compare && i < n; i += 1)

            // add a[i] to the elements that have been gap sorted
            // save a[i] in temp and make a hole at position i
        {
            int temp = arr[i];
            // shift earlier gap-sorted elements up until the correct 
           // location for a[i] is found
            int j;
            for (j = i; ++cnt_compare && j >= gap && ++cnt_compare && arr[j - gap] > temp; j -= gap)
                arr[j] = arr[j - gap];
            //  put temp (the original a[i]) in its correct location
            arr[j] = temp;
        }
    }
}

int getMax(int arr[], int n)
{
    int mx = arr[0];
    for (int i = 1; ++cnt_compare && i < n; i++)
        if (++cnt_compare && arr[i] > mx)
            mx = arr[i];
    return mx;
}
void countSort(int arr[], int n, int exp)
{
    int* output = new int[n];
    int i, count[10] = { 0 };

    for (i = 0; ++cnt_compare && i < n; i++)
        count[(arr[i] / exp) % 10]++;

    for (i = 1; ++cnt_compare && i < 10; i++)
        count[i] += count[i - 1];

    for (i = n - 1; ++cnt_compare && i >= 0; i--) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
    }
    for (i = 0; ++cnt_compare && i < n; i++)
        arr[i] = output[i];
}
void radixSort(int arr[], int n)
{
    int m = getMax(arr, n);
    for (int exp = 1; ++cnt_compare && m / exp > 0; exp *= 10)
        countSort(arr, n, exp);
}


void selectionSort(int a[], int n)// tham khảo từ github
{

    for (int i = 0; ++cnt_compare, i < n - 1; i++)

    {
        int vitriMin = i;

        for (int j = i + 1; ++cnt_compare, j < n; j++)
        {
            if (++cnt_compare && a[j] < a[vitriMin])
            {
                vitriMin = j;
            }
        }
        HoanVi(a[i], a[vitriMin]);
    }
}

void merge(int a[], int left, int mid, int right)// tham khảo từ geeksforgeeks
{
    int* b = new int[right - left + 1]; // m?ng t?m l?u m?ng con sau khi x?p
    int c1 = 0, c2 = 0;
    int d = 0;

    for (d = 0, c1 = left, c2 = mid + 1; (c1 <= mid) && (c2 <= right); d++) // x�t c�c ph?n t? t? left->mid v� mid+1->right
    {
        if (a[c1] > a[c2])
        {
            b[d] = a[c1];
            c1++;
        }
        else
        {
            b[d] = a[c2];
            c2++;
        }
    }

    while (c1 <= mid) // n?u m?ng t? left -> mid c�n ph?n t? ta g�n l?n l??t v�o m?ng kqua.
    {
        b[d++] = a[c1++];
    }

    while (c2 <= right) // n?u m?ng t? mid+1 -> right c�n ph?n t? ta g�n l?n l??t v�o m?ng kqua
    {
        b[d++] = a[c2++];
    }

    for (int i = 0; i < d; i++) //g�n l?i c�c ph?n t? c?a m?ng k?t qu? v�o m?ng ban ??u
    {
        a[i + left] = b[i];
    }
}

void mergeSort(int arr[], int l, int r) {
    if (++cnt_compare && l < r) {
        // m is the point where the array is divided into two subarrays
        int m = l + (r - l) / 2;

        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        // Merge the sorted subarrays
        merge(arr, l, m, r);
    }
}

void mergeSort(int a[], int n)
{
    mergeSort(a, 0, n - 1);
}

void flashSort(int a[], int n) // tham khảo từ github
{
    int minVal = a[0];
    int max = 0;
    int m = int(0.45 * n);
    int* l = new int[m];
    for (int i = 0; ++cnt_compare, i < m; i++)
        l[i] = 0;
    for (int i = 1; ++cnt_compare, i < n; i++)
    {
        if (++cnt_compare && a[i] < minVal)
            minVal = a[i];
        if (++cnt_compare && a[i] > a[max])
            max = i;
    }
    if (++cnt_compare && a[max] == minVal)
        return;
    double c1 = (double)(m - 1) / (a[max] - minVal);
    for (int i = 0; ++cnt_compare, i < n; i++)
    {
        int k = int(c1 * (a[i] - minVal));
        ++l[k];
    }
    for (int i = 1; ++cnt_compare, i < m; i++)
        l[i] += l[i - 1];
    HoanVi(a[max], a[0]);
    int nmove = 0;
    int j = 0;
    int k = m - 1;
    int t = 0;
    int flash;
    while (nmove < n - 1)
    {
        while (j > l[k] - 1)
        {
            j++;
            k = int(c1 * (a[j] - minVal));
        }
        flash = a[j];
        if (k < 0) break;
        while (j != l[k])
        {
            k = int(c1 * (flash - minVal));
            int hold = a[t = --l[k]];
            a[t] = flash;
            flash = hold;
            ++nmove;
        }
    }
    delete[] l;
    insertionSort(a, n);
}

void printTimeSpent(double time_spent)
{
    //cout << "Running time: " << time_spent << " s" << endl;
    printf("Running time: %0.3lf mms\n", time_spent*1000);
}
void printComparison(unsigned long long int compare)
{
    cout << "Comparisons: " << compare << endl;
}
void printTimeSpent(clock_t start, clock_t end)
{
    long double time_use = (long double)(end - start) / CLOCKS_PER_SEC;
    cout << "Running time: " << time_use << " second.\n";
}
void cout_random()
{
    cout << "\nInput order: Randomize" << endl;
    cout << "-------------------------" << endl;
}
void cout_sorted()
{
    cout << "\nInput order: Sorted" << endl;
    cout << "-------------------------" << endl;
}
void cout_nearlysorted()
{
    cout << "\nInput order: Nearly Sorted" << endl;
    cout << "-------------------------" << endl;
}
void cout_reversed()
{
    cout << "\nInput order: Reversed" << endl;
    cout << "-------------------------" << endl;
}

void (*function_algorithm[11])(int*, int) = { &selectionSort,&insertionSort,&bubbleSort,&shakerSort,&shellSort,&heapSort,&mergeSort,&quickSort,&countingSort,&radixSort,&flashSort };
string algorithm_name[11] = { "selection-sort", "binary-insertion-sort", "bubble-sort","shaker-sort","shell-sort","heap-sort","merge-sort","quick-sort","counting-sort","radix-sort","flash-sort" };

void calculate_algorithm(string algorithm, clock_t begin, double& time_spent)
{
    for (int i = 0; i < 11; i++)
    {
        if (algorithm_name[i] == algorithm)
        {
            function_algorithm[i](a, n);
        }
    }

    clock_t end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
}

void command_1(string action, string algorithm, string input_file, string output_parameter)
{
    ifstream fin;
    fin.open(input_file);
    fin >> n;
    for (int i = 0; i < n; i++)
        fin >> a[i];
    fin.close();

    cout << "ALGORITHM MODE" << endl;
    cout << "Algorithm: " << algorithm << endl;
    cout << "Input file: " << input_file << endl;
    cout << "Input size: " << n << endl;

    clock_t begin = clock();
    double time_spent = 0;
    calculate_algorithm(algorithm, begin, time_spent);
    //cout << begin << ' ' << end << endl;
//printf("%0.3lf\n", time_spent);
    cout << "-------------------------" << endl;
    if (output_parameter == "-time")
    {
        printTimeSpent(time_spent);
    }
    else if (output_parameter == "-comp")
    {
        printComparison(cnt_compare);
    }
    else
    {
        printTimeSpent(time_spent);
        printComparison(cnt_compare);
    }

    ofstream fout;
    fout.open("output.txt");
    fout << n << endl;
    for (int i = 0; i < n; i++)
        fout << a[i] << ' ';
    fout.close();
}
void command_2(string algorithm, int input_size, string input_order, string output_parameter)
{
    //cout << algorithm << ' ' << input_size << ' ' << input_order << ' ' << output_parameter << endl;
    n = input_size;
    if (input_order == "-rand")
        GenerateData(a, n, 0);
    if (input_order == "-nsorted")
        GenerateData(a, n, 3);
    if (input_order == "-sorted")
        GenerateData(a, n, 1);
    if (input_order == "-rev")
        GenerateData(a, n, 2);

    ofstream fout;
    fout.open("input.txt");
    fout << n << endl;
    for (int i = 0; i < n; i++)
        fout << a[i] << ' ';
    fout.close();

    clock_t begin = clock();
    cout << endl << "ALGORITHM MODE" << endl;
    cout << "Algorithm: " << algorithm << endl;
    cout << "Input size: " << n << endl;
    cout << "Input order: " << input_order << endl;

    double time_spent = 0;
    calculate_algorithm(algorithm, begin, time_spent);

    cout << "-------------------------" << endl;
    if (output_parameter == "-time")
    {
        printTimeSpent(time_spent);
    }
    else if (output_parameter == "-comp")
    {
        printComparison(cnt_compare);
    }
    else
    {
        printTimeSpent(time_spent);
        printComparison(cnt_compare);
    }

    fout.open("output.txt");
    fout << n << endl;
    for (int i = 0; i < n; i++)
        fout << a[i] << ' ';
    fout.close();
}
void command_3(string mode, string algorithm, string input_size, string output_parameter)
{
    if (mode == "-a") cout << "ALGORITHM MODE\n";
    else cout << "COMPARISION MODE\n";
    int size;
    clock_t start, end;
    int* a = new int[500000];
    size = stoi(input_size);

    if (algorithm == "heap-sort")
    {
        GenerateRandomData(a, size);
        cout_random();
        write_input(a, size, 1);
        start = clock();
        heapSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateNearlySortedData(a, size);
        cout_nearlysorted();
        write_input(a, size, 2);
        start = clock();
        heapSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateSortedData(a, size);
        cout_sorted();
        write_input(a, size, 3);
        start = clock();
        heapSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateReverseData(a, size);
        cout_reversed();
        write_input(a, size, 4);
        start = clock();
        heapSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
    }
    else if (algorithm == "binary-insertion-sort")
    {
        GenerateRandomData(a, size);
        cout_random();
        write_input(a, size, 1);
        start = clock();
        insertionSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateNearlySortedData(a, size);
        cout_nearlysorted();
        write_input(a, size, 2);
        start = clock();
        insertionSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateSortedData(a, size);
        cout_sorted();
        write_input(a, size, 3);
        start = clock();
        insertionSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateReverseData(a, size);
        cout_reversed();
        write_input(a, size, 4);
        start = clock();
        insertionSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
    }
    else if (algorithm == "bubble-sort")
    {
        GenerateRandomData(a, size);
        cout_random();
        write_input(a, size, 1);
        start = clock();
        bubbleSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateNearlySortedData(a, size);
        cout_nearlysorted();
        write_input(a, size, 2);
        start = clock();
        bubbleSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateSortedData(a, size);
        cout_sorted();
        write_input(a, size, 3);
        start = clock();
        bubbleSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateReverseData(a, size);
        cout_reversed();
        write_input(a, size, 4);
        start = clock();
        bubbleSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
    }
    else if (algorithm == "quick-sort")
    {
        GenerateRandomData(a, size);
        cout_random();
        write_input(a, size, 1);
        cnt_compare = 0;
        start = clock();
        QuickSort(a, 0, size - 1);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateNearlySortedData(a, size);
        cout_nearlysorted();
        write_input(a, size, 2);
        cnt_compare = 0;
        start = clock();
        QuickSort(a, 0, size - 1);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateSortedData(a, size);
        cout_sorted();
        write_input(a, size, 3);
        cnt_compare = 0;
        start = clock();
        QuickSort(a, 0, size - 1);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateReverseData(a, size);
        cout_reversed();
        write_input(a, size, 4);
        cnt_compare = 0;
        start = clock();
        QuickSort(a, 0, size - 1);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
    }
    else if (algorithm == "shell-sort")
    {
        GenerateRandomData(a, size);
        cout_random();
        write_input(a, size, 1);
        start = clock();
        shellSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateNearlySortedData(a, size);
        cout_nearlysorted();
        write_input(a, size, 2);
        start = clock();
        shellSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateSortedData(a, size);
        cout_sorted();
        write_input(a, size, 3);
        start = clock();
        shellSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateReverseData(a, size);
        cout_reversed();
        write_input(a, size, 4);
        start = clock();
        shellSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
    }
    else if (algorithm == "selection-sort")
    {
        GenerateRandomData(a, size);
        cout_random();
        write_input(a, size, 1);
        start = clock();
        selectionSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateNearlySortedData(a, size);
        cout_nearlysorted();
        write_input(a, size, 2);
        start = clock();
        selectionSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateSortedData(a, size);
        cout_sorted();
        write_input(a, size, 3);
        start = clock();
        selectionSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateReverseData(a, size);
        cout_reversed();
        write_input(a, size, 4);
        start = clock();
        selectionSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
    }
    else if (algorithm == "merge-sort")
    {
        GenerateRandomData(a, size);
        cout_random();
        write_input(a, size, 1);
        cnt_compare = 0;
        start = clock();
        mergeSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateNearlySortedData(a, size);
        cout_nearlysorted();
        write_input(a, size, 2);
        cnt_compare = 0;
        start = clock();
        mergeSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateSortedData(a, size);
        cout_sorted();
        write_input(a, size, 3);
        cnt_compare = 0;
        start = clock();
        mergeSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateReverseData(a, size);
        cout_reversed();
        write_input(a, size, 4);
        cnt_compare = 0;
        start = clock();
        mergeSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
    }
    else if (algorithm == "shaker-sort")
    {
        GenerateRandomData(a, size);
        cout_random();
        write_input(a, size, 1);
        start = clock();
        shakerSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateNearlySortedData(a, size);
        cout_nearlysorted();
        write_input(a, size, 2);
        start = clock();
        shakerSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateSortedData(a, size);
        cout_sorted();
        write_input(a, size, 3);
        start = clock();
        shakerSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateReverseData(a, size);
        cout_reversed();
        write_input(a, size, 4);
        start = clock();
        shakerSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
    }
    else if (algorithm == "radix-sort")
    {
        GenerateRandomData(a, size);
        cout_random();
        write_input(a, size, 1);
        start = clock();
        radixSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateNearlySortedData(a, size);
        cout_nearlysorted();
        write_input(a, size, 2);
        start = clock();
        radixSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateSortedData(a, size);
        cout_sorted();
        write_input(a, size, 3);
        start = clock();
        radixSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
        GenerateReverseData(a, size);
        cout_reversed();
        write_input(a, size, 4);
        start = clock();
        radixSort(a, size);
        end = clock();
        if (output_parameter == "-time" || output_parameter == "-both") printTimeSpent(start, end);
        if (output_parameter == "-comp" || output_parameter == "-both") printComparison(cnt_compare);
    }
}
void command_4(string action, string algorithm1, string algorithm2, string input_file)
{
    ifstream fin;
    fin.open(input_file);
    fin >> n;
    int* temparray = new int[n];
    for (int i = 0; i < n; i++)
    {
        fin >> a[i];
        temparray[i] = a[i];
    }
    fin.close();

    clock_t begin = clock();
    double time_spent1 = 0;

    calculate_algorithm(algorithm1, begin, time_spent1);
    unsigned long long int compare1 = cnt_compare;

    clock_t begin2 = clock();
    double time_spent2 = 0;
    cnt_compare = 0;
    calculate_algorithm(algorithm2, begin2, time_spent2);
    unsigned long long int compare2 = cnt_compare;

    cout << endl << "COMPARE MODE" << endl;
    cout << "Algorithm: " << algorithm1 << " | " << algorithm2 << endl;
    cout << "Input file: " << input_file << endl;
    cout << "Input size: " << n << endl;
    cout << "-------------------------" << endl;

    cout << "Running time: " << time_spent1 << " | " << time_spent2 << endl;
    cout << "Comparisions: " << compare1 << " | " << compare2 << endl;
}
void command_5(string algorithm1, string algorithm2, int input_size, string input_order)
{
    n = input_size;

    if (input_order == "-rand")
        GenerateData(a, n, 0);
    if (input_order == "-nsorted")
        GenerateData(a, n, 3);
    if (input_order == "-sorted")
        GenerateData(a, n, 1);
    if (input_order == "-rev")
        GenerateData(a, n, 2);
    /*   for (int i = 0; i < n; i++)
           cout << a[i] << ' ';
       cout << endl;*/
    ofstream fout;
    fout.open("input.txt");
    fout << n << endl;
    for (int i = 0; i < n; i++)
        fout << a[i] << ' ';
    fout.close();

    clock_t begin = clock();
    double time_spent1 = 0;

    calculate_algorithm(algorithm1, begin, time_spent1);
    unsigned long long int compare1 = cnt_compare;

    clock_t begin2 = clock();
    double time_spent2 = 0;
    cnt_compare = 0;
    calculate_algorithm(algorithm2, begin2, time_spent2);
    unsigned long long int compare2 = cnt_compare;

    cout << endl << "COMPARE MODE" << endl;
    cout << "Algorithm: " << algorithm1 << " | " << algorithm2 << endl;
    cout << "Input size: " << n << endl;
    cout << "Input order: " << input_order << endl;
    cout << "-------------------------" << endl;

    cout << "Running time: " << time_spent1 << " | " << time_spent2 << endl;
    cout << "Comparisions: " << compare1 << " | " << compare2 << endl;
}

int main(int argc, char* argv[])
{
    string input_file, action, algorithm, output_parameter, input_order;
    int input_size = 0;
    
   if (argv[1][1] == 'a')
   {
        if (argc == 6)
        {
            action = argv[1];
            algorithm = argv[2];
            input_size = stoi(argv[3]);
            input_order = argv[4];
            output_parameter = argv[5];
            //cout << "cmd2" << endl;
            command_2(algorithm, input_size, input_order, output_parameter);
            return 0;
        }
        else
        {
            if (isNum(argv[3]))
            {
                action = argv[1];
                algorithm = argv[2];
                string input_size = argv[3];
                output_parameter = argv[4];
                //cout << "cmd3" << endl;
                command_3(action, algorithm, input_size, output_parameter);
                return 0;
            }
            else
            {
                action = argv[1];
                algorithm = argv[2];
                input_file = argv[3];
                output_parameter = argv[4];
                //cout << "cmd1" << endl;
                command_1(action, algorithm, input_file, output_parameter);
                return 0;
            }
        }
    }
    else
    {
        if (argc == 5)
        {
            action = argv[1];
            string algorithm1 = argv[2];
            string algorithm2 = argv[3];
            input_file = argv[4];
            //cout << "cmd4" << endl;
            command_4(action, algorithm1, algorithm2, input_file);
            return 0;
        }
        else
        {
            action = argv[1];
            string algorithm1 = argv[2];
            string algorithm2 = argv[3];
            input_size = stoi(argv[4]);
            input_order = argv[5];
            //cout << "cmd5" << endl;
            command_5(algorithm1, algorithm2, input_size, input_order);
            return 0;
        }
    }
    if (!system(NULL)) system("pause"); return 0;
    return 0;
}
