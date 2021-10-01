#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>

using namespace std;

struct PLAYER 
{
    string name;
    int money;

    PLAYER(string name)
    {
        this->name = name;
        money = 0;
    }
};

unordered_map<string, PLAYER> players;

void enters(string& name)
{
    unordered_map<string, PLAYER>::iterator player = players.find(name);
    if (player != players.end()) 
    {
        cout << name << " is already in the casino" << endl;
        return;
    }

    players.insert(std::pair<string, PLAYER>(name, PLAYER(name)));

}

void wins(string& name)
{
    unordered_map<string, PLAYER>::iterator player = players.find(name);
    int coins; cin >> coins;
    if (player == players.end())
    {
        cout << name << " is not in the casino" << endl;
        return;
    }

    player->second.money += coins;
}

bool cmp (unordered_map<string, PLAYER>::const_iterator& a, unordered_map<string, PLAYER>::const_iterator& b)
{
    return (*a).second.name < (*b).second.name;
}

void leaves(string& name)
{
    unordered_map<string, PLAYER>::iterator player = players.find(name);
    if (player == players.end())
    {
        cout << name << " is not in the casino" << endl;
        return;
    }

    cout << player->second.name << " has won " << player->second.money << endl;
    players.erase(player);
}

int main ()
{
    string name;
    while (cin >> name)
    {
        string action; cin >> action;

        if (action == "enters"  ) enters    (name);
        if (action == "wins"    ) wins      (name);
        if (action == "leaves"  ) leaves    (name);

    }
    cout << "----------" << endl;
    vector<unordered_map<string, PLAYER>::const_iterator> vect(players.size());
    unordered_map<string, PLAYER>::iterator iter = players.begin();

    for (int i = 0; i < vect.size(); ++i) vect[i] = iter, ++iter;

    sort(vect.begin(), vect.end(), cmp);

    for (int i = 0; i < vect.size(); ++i)
        cout << vect[i]->second.name << " is winning " << vect[i]->second.money  << endl;
}