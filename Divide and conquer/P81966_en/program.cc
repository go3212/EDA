#include <iostream>
#include <vector>

using namespace std;

int position(double x, const vector<double>& v, int left, int right)
{
    // Basicamente queremos que el algoritmo se ejecute de la siguiente manera
    // 1. 0 = L, R = (size - 1) *nos llega por position*
    // 2. Si L > R acabar busqueda *not found* return -1;
    // 3. Sea m el elemento medio entre L y R... m = (L + R) / 2
    // 4. Si v[m] < x -> L = m + 1 ---> paso 2.
    // 5. Si v[m] > x -> R = m - 1 ---> paso 2
    // 6. Si v[m] == x -> return x;

    if (left > right) return -1;
    int m = (left + right) / 2;
    if (v[m] < x) return position(x, v, m + 1, right);
    if (v[m] > x) return position(x, v, left, m - 1);
    return m;

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
            int left, right;
            cin >> x >> left >> right;
            cout << position(x, V, left, right) << endl;
        }
    }
}
