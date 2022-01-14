#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

// S = x_1, ..., x_n tal que x_1 < ... < x_n. 

int getFixed(int& a, vector<int>& S)
{
    // Implemento un hibrido de búsqueda dicotómica. x_i + a < x_k + a, si i < k (por definición).
    int left = 0, right = S.size() - 1, mid;
    while (left <= right)
    {
        mid = (left+right)/2; 
        int i = S[mid] + a;

        // Primero eliminamos la posibilidad
        if (i > mid + 1)                                        {right = mid - 1; continue;};
        if (i < mid + 1)                                        {left  = mid + 1; continue;};
        if (i == mid + 1 && mid > 0 && S[mid - 1] + a == mid)   {right = mid - 1; continue;};
        return mid+1;
    }
    return -1;
}

int main ()
{
	int n, seq;
	seq = 1;
	while (cin >> n) {
		cout << "Sequence #" << seq << endl;
		vector<int> v(n);
		int m;
		for(int i = 0; i < n; ++i) cin >> v[i];
		cin >> m;
		for (int i = 0; i < m; ++i) {
			int a, fixed;
			cin >> a;
			fixed = getFixed(a, v);
			if (fixed == -1) cout << "no fixed point for " << a;
			else cout << "fixed point for " << a << ": " << fixed;
			cout << endl;
		}
		++seq;
		cout << endl;
	}
}