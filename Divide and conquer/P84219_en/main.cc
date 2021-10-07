#include <iostream>
#include <vector>

using namespace std;

int first_occurrence(double x, const vector<double>& v)
{
    if (v.size() == 0) return -1;
    // El vector está ordenado en orden creciente. V[i - 1] < V[i]. Empleamos búsqueda binária para el primer elemento.
    int left = 0, right = v.size() - 1;
    int mid = 0;
    while (left <= right)
    {
        mid = (left+right)/2;
        if (v[mid] < x) { left  = mid + 1; continue; }
        if (v[mid] > x) { right = mid - 1; continue; }
        if (mid > 0 && v[mid - 1] == x) {right = mid - 1; continue;}
        return mid;
    }
    return -1;
}


int main() {
    int n;
    while (cin >> n) {
        vector<double> V(n);
        for (int i = 0; i < n; ++i) cin >> V[i];
        int t;
        cin >> t;
        while (t--) {
            double x;
            cin >> x;
            cout << first_occurrence(x, V) << endl;
        }
    }
}
// Casos:
// 6 1 2 2 2 3 3 1 2 -- > 1