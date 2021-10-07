#include <iostream>

using namespace std;

int mod_exp(int n, int k, int m)
{
    if (k == 0) return 1;
    if (k % 2 != 0) return (n*mod_exp(n, k-1, m))%m;
    int temp = mod_exp(n, k/2, m)%m;
    return (temp*temp)%m;

}

int main ()
{
    // Programa que dado n, k y m compute (n^k mod m)
    // n^k % m --> n*n*n*n*...*n_k
    int n, k, m;
    while (cin >> n >> k >> m)
    {
        cout << mod_exp(n, k, m) << endl;
    }   
}