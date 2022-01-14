#include <iostream>
#include <vector>

int* max_value (std::vector<int>& vect)
{
    int size = vect.size();
    if (size == 0) return NULL;
    
    int* max = &vect[0];
    for (int i = 1; i < size; ++i)
        max = ((vect[i] >= *max) ? &vect[i] : max);

    return max;
}

int main ()
{
    int* max = NULL;
    std::vector<int> storage;

    char input;
    while (std::cin >> input)
    {
        int x; 
        switch (input)
        {
            case 'S':
                std::cin >> x;
                storage.push_back(x);
                if (max == NULL) max = &storage[storage.size() - 1];
                else max = ((x >= *max) ? &storage[storage.size() - 1] : max);
                
                break;
            case 'A':
                std::cout << ((max == NULL) ? ("error!") : std::to_string(*max)) << std::endl; 
                break;
            case 'R':
                if (max == NULL) 
                {
                    std::cout << "error!" << std::endl;
                    break;
                }
                // *max es el valor mas grande, recalculamos la posicion
                storage.erase(storage.begin() + (max - &storage[0])); // funciona ya que al menos tenemos un elemento.
                // Ahora buscamos el valor máximo.        
                max = max_value(storage);
                break;
            case 'I':
                if (max == NULL) 
                {
                    std::cout << "error!" << std::endl;
                    break;
                }
                std::cin >> x;
                *max += x;
                max = max_value(storage);
                break;
            case 'D':
                if (max == NULL) 
                {
                    std::cout << "error!" << std::endl;
                    break;
                }
                std::cin >> x;
                *max -= x;
                max = max_value(storage);
                break;

        }

    }
}