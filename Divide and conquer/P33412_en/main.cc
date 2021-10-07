#include <iostream>
#include <vector>
using namespace std;


bool resistant_search(double x, const vector<double>& v)
{
    int left, right = v.size() - 1, mid;
    while(left <= right)
    {
        mid = (left+right)/2;
        if (v[mid] == x) return true;
        if (mid < v.size() - 1 && v[mid] > v[mid + 1] && v[mid + 1] == x) return true;
        if (mid > 0            && v[mid] < v[mid - 1] && v[mid - 1] == x) return true;
        if (v[mid] > x) { right = mid - 1; continue;}
        if (v[mid] < x) { left  = mid + 1; continue;}
    }
    return false;
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
            cout << resistant_search(x, V) << endl;
        }
    }
}
// 6 0 1 2 3 4 5 1 2 --> 1
// 6 1 0 2 3 4 5 1 1 --> 1
// 6 0 1 2 3 5 4 1 4 --> 1
// 6 0 1 2 3 5 4 1 5 --> 1
// 6 0 1 3 2 4 5 1 5 --> 1