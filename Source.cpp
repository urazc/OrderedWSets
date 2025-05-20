#include<fstream>
#include<iostream>
#include<vector>
#include <stack>
#include <queue>
#include <chrono>
#include <Windows.h>
#include <Psapi.h>
#include <deque>
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
using namespace std;
vector<vector<int>> nextStates;
vector<vector<int>> observedOutputs;
vector<vector<vector<int>>> separatingSequences;
int FSMid = 0, FSMStates = 0, FSMSize = 0, FSMInputs = 0, FSMOutputs = 0, dv = 0;
int pairCount = 0;
int transferLengthC = 0;
int transferLengthO = 0;
int transferLengthM = 0;

int maximumPath = 0;

int depth = 0;
deque<int> indexes;
vector<int> inputsOrder;
stack<int> mySequence;
stack<int> sequence;
struct edge {
    int dest;
    bool visited;
};

vector<vector<edge>> edges;
queue<pair<int, vector<int>>> myQ;
vector<vector<bool>> sms;


double getUsedMemoryMB()//Get the peak working set size (i.e., peak memory usage) of the current process in MB.
{
    PROCESS_MEMORY_COUNTERS_EX pmc;
    DWORD ret = GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
    if (ret == 0)
    {
        printf("GetProcessMemoryInfo failed with code %d\n", GetLastError());
    }
    size_t temp2 = pmc.PeakWorkingSetSize;
    double mem = static_cast<float>(temp2);
    mem = mem / 1048576.0;

    return mem;
}
unsigned int returnPairIndex(int i, int j, int n)//Compute a unique index for a pair (i, j) where i < j and n is the total number of elements.
{
    int result = ((i * n) - (i * (i + 1) / 2)) + (j - i - 1);
    return result;
}
vector<pair<int, int>> separates(vector<int>& s, int x)//Separate states based on their output for input x and return pairs of states that can be separated.
{
    vector<vector<int>>temp;
    temp.resize(FSMOutputs);

    for (int i = 0; i < s.size(); i++)
    {
        temp[observedOutputs[s[i]][x]].push_back(s[i]);
    }
    vector<vector<int>>states;
    for (int i = 0; i < temp.size(); i++)
    {
        if (temp[i].size() == 0)
            ;
        else
        {
            states.push_back(temp[i]);
        }
    }

    vector<pair<int, int>> origins;
    if (states.size() == 1)
        return origins;
    vector<int> seq;
    seq.push_back(x);
    for (int i = 0; i < states.size() - 1; i++)
    {
        for (int j = i + 1; j < states.size(); j++)
        {
            for (int k = 0; k < states[i].size(); k++)
            {
                for (int l = 0; l < states[j].size(); l++)
                {
                    if (states[i][k] < states[j][l])
                    {
                        origins.push_back(pair<int, int>(states[i][k], states[j][l]));

                        if (separatingSequences[returnPairIndex(states[i][k], states[j][l], observedOutputs.size())].size() == 0)
                            pairCount--;
                        separatingSequences[returnPairIndex(states[i][k], states[j][l], observedOutputs.size())].push_back(seq);

                    }
                    else
                    {

                        origins.push_back(pair<int, int>(states[j][l], states[i][k]));
                        if (separatingSequences[returnPairIndex(states[j][l], states[i][k], observedOutputs.size())].size() == 0)
                            pairCount--;
                        separatingSequences[returnPairIndex(states[j][l], states[i][k], observedOutputs.size())].push_back(seq);
                    }
                }
            }
        }
    }
    return origins;
}
vector<pair<int, int>> separatesClassic(vector<int>& s, int x)//Classic version of separates, with slightly different sequence handling logic.
{
    vector<vector<int>>temp;
    temp.resize(FSMOutputs);

    for (int i = 0; i < s.size(); i++)
    {
        temp[observedOutputs[s[i]][x]].push_back(s[i]);
    }
    vector<vector<int>>states;
    for (int i = 0; i < temp.size(); i++)
    {
        if (temp[i].size() == 0)
            ;
        else
        {
            states.push_back(temp[i]);
        }
    }

    vector<pair<int, int>> origins;
    if (states.size() == 1)
        return origins;
    vector<int> seq;
    seq.push_back(x);
    for (int i = 0; i < states.size() - 1; i++)
    {
        for (int j = i + 1; j < states.size(); j++)
        {
            for (int k = 0; k < states[i].size(); k++)
            {
                for (int l = 0; l < states[j].size(); l++)
                {
                    if (states[i][k] < states[j][l])
                    {
                        origins.push_back(pair<int, int>(states[i][k], states[j][l]));

                        if (separatingSequences[returnPairIndex(states[i][k], states[j][l], observedOutputs.size())].size() == 0) {
                            pairCount--;
                            separatingSequences[returnPairIndex(states[i][k], states[j][l], observedOutputs.size())].push_back(seq);
                        }

                    }
                    else
                    {

                        origins.push_back(pair<int, int>(states[j][l], states[i][k]));
                        if (separatingSequences[returnPairIndex(states[j][l], states[i][k], observedOutputs.size())].size() == 0) {
                            pairCount--;
                            separatingSequences[returnPairIndex(states[j][l], states[i][k], observedOutputs.size())].push_back(seq);
                        }
                    }
                }
            }
        }
    }
    return origins;
}
pair< vector<pair<int, int>>, vector<int>> incrementAClass(pair< vector<pair<int, int>>, vector<int>>& myClass, int input)//Generate new distinguishing pairs from an existing class and an input, updating the separating sequences.
{
    vector<pair<int, int>> RESULT;
    vector<int> seq = myClass.second;
    seq.push_back(input);
    //seq.push_back(input);

    if (myClass.second.size() == 0)
    {
        return  pair< vector<pair<int, int>>, vector<int>>(vector<pair<int, int>>(), vector<int>());
    }
    /*for (int i = 0; i < myClass.second.size(); i++)
    {
        seq.push_back(myClass.second[i]);
    }*/

    for (int i = 0; i < myClass.first.size(); i++)
    {
        int L = myClass.first[i].first;
        int R = myClass.first[i].second;


        vector<int> tempL, tempR;

        for (int k = 0; k < nextStates.size(); k++)
        {
            if (L != R)
            {
                {
                    if (nextStates[k][input] == L)
                    {
                        tempL.push_back(k);

                    }
                    if (nextStates[k][input] == R)
                    {
                        tempR.push_back(k);
                    }
                }
            }

        }

        if (tempR.size() > 0 && tempL.size() > 0)
        {

            for (int k = 0; k < tempL.size(); k++)
            {
                for (int l = 0; l < tempR.size(); l++)
                {
                    if (tempL[k] >= 0 && tempR[l] >= 0 && tempL[k] < tempR[l])
                    {
                        RESULT.push_back(pair<int, int>(tempL[k], tempR[l]));
                        int index = returnPairIndex(tempL[k], tempR[l], observedOutputs.size());

                        if (separatingSequences[index].size() == 0)
                        {
                            pairCount--;
                            separatingSequences[index].push_back(seq);
                        }
                        else 
                        {
                            //separatingSequences[index].push_back(seq);
                        }

                    }
                    else if (tempL[k] >= 0 && tempR[l] >= 0 && tempL[k] > tempR[l])
                    {
                        RESULT.push_back(pair<int, int>(tempR[l], tempL[k]));
                        int index = returnPairIndex(tempR[l], tempL[k], observedOutputs.size());

                        if (separatingSequences[index].size() == 0)
                        {
                            pairCount--;
                            separatingSequences[index].push_back(seq);
                        }
                        else //if (input != seq[0])
                        {
                            //separatingSequences[index].push_back(seq);
                        }

                    }
                }
            }

        }
    }
    if (RESULT.size() > 0)
        return pair< vector<pair<int, int>>, vector<int>>(RESULT, seq);
    else
        return  pair< vector<pair<int, int>>, vector<int>>(vector<pair<int, int>>(), vector<int>());
}
pair< vector<pair<int, int>>, vector<int>> incrementAClassClassic(pair< vector<pair<int, int>>, vector<int>>& myClass, int input)//Classic version of incrementAClass with slightly different logic.
{
    vector<pair<int, int>> RESULT;
    vector<int> seq;
    seq.push_back(input);

    for (int i = 0; i < myClass.second.size(); i++)
    {
        seq.push_back(myClass.second[i]);
    }

    for (int i = 0; i < myClass.first.size(); i++)
    {
        int L = myClass.first[i].first;
        int R = myClass.first[i].second;

        if (L == 1 && R == 4 && input == 1)
            int axd = 0;
        vector<int> tempL, tempR;

        for (int k = 0; k < nextStates.size(); k++)
        {
            if (L != R)
            {
                {
                    if (nextStates[k][input] == L)
                    {
                        tempL.push_back(k);
                    }
                    if (nextStates[k][input] == R)
                    {
                        tempR.push_back(k);
                    }
                }
            }

        }

        if (tempR.size() > 0 && tempL.size() > 0)
        {

            for (int k = 0; k < tempL.size(); k++)
            {
                for (int l = 0; l < tempR.size(); l++)
                {
                    if (tempL[k] >= 0 && tempR[l] >= 0 && tempL[k] < tempR[l])
                    {
                        RESULT.push_back(pair<int, int>(tempL[k], tempR[l]));
                        int index = returnPairIndex(tempL[k], tempR[l], observedOutputs.size());

                        if (separatingSequences[index].size() == 0)
                        {
                            pairCount--;
                            separatingSequences[index].push_back(seq);
                        }


                    }
                    else if (tempL[k] >= 0 && tempR[l] >= 0 && tempL[k] > tempR[l])
                    {
                        RESULT.push_back(pair<int, int>(tempR[l], tempL[k]));
                        int index = returnPairIndex(tempR[l], tempL[k], observedOutputs.size());

                        if (separatingSequences[index].size() == 0)
                        {
                            pairCount--;
                            separatingSequences[index].push_back(seq);
                        }


                    }
                }
            }

        }
    }
    if (RESULT.size() > 0)
        return pair< vector<pair<int, int>>, vector<int>>(RESULT, seq);
    else
        return  pair< vector<pair<int, int>>, vector<int>>(vector<pair<int, int>>(), vector<int>());
}
bool isAdded(vector<int>& l, int r)//Check if value r is already in list l.
{
    bool nt = false;
    for (int i = 0; i < l.size(); i++)
    {
        if (l[i] == r)
            return true;

    }
    return false;
}
bool isPrefix(vector<int>& L, vector<int>& R)//Check if vector L is a prefix of vector R.
{
    for (int i = 0; i < L.size(); i++)
    {
        if (L[i] != R[i])
        {
            return false;
        }

    }
    return true;
}
bool areEqual(vector<int>& L, vector<int>& R)//Check if two sequences are equal up to the shorter one being a prefix.
{
    if (L.size() < R.size())
    {
        if (isPrefix(L, R))
            return true;
    }
    else if (R.size() < L.size())
    {
        if (isPrefix(R, L))
            return true;
    }
    return isPrefix(L, R);
}
bool isLinR(vector<int>& L, vector<int>& R)//Check if L is a prefix (or equal to) of R.
{
    if (L.size() <= R.size())
    {
        if (isPrefix(L, R))
            return true;
    }

    return false;
}
int isIn(vector<int>& L, vector<int>& R)//Check if L is a prefix of R and return 1, else 0, or -1 if L is longer than R.
{
    if (L.size() > R.size())
    {
        return -1;
    }
    else
    {
        for (int i = 0; i < L.size(); i++)
        {
            if (L[i] != R[i])
                return 0;
        }
    }
    return 1;
}
bool hasAsset(vector<vector<int>>& cs, vector<int>& seqs)//Check if a given sequence seqs exists in the condensed set cs.
{
    for (int i = 0; i < cs.size(); i++)
    {
        int r = isIn(seqs, cs[i]);
        if (r == 1)
            return true;
    }
    return false;
}
int reachedState(int state, vector<int> seq)//Determine the final state reached from a start state following a sequence.
{
    for (int i = 0; i < seq.size(); i++)
    {
        state = nextStates[state][seq[i]];
    }
    return state;
}
int getOrderingSequence(vector<vector<int>>& seqs)//Select the sequence that best distinguishes all states, prioritizing minimal length.
{
    int max_transfer_index = INT_MIN;
    int traLength = INT_MAX;
    vector<int> transferIndexes;
    vector<int> transferIndexesSize;

    for (int i = 0; i < seqs.size(); i++)
    {
        vector<int> Targets;
        int transferIndex = 0;
        Targets.resize(nextStates.size());
        Targets.assign(nextStates.size(), 0);
        for (int j = 0; j < nextStates.size(); j++)
        {
            int icIndex = reachedState(j, seqs[i]);
            Targets[icIndex]++;
        }
        int f = 0;
        for (int j = 0; j < nextStates.size(); j++)
        {
            if (Targets[j] == 1)
                transferIndex++;
        }
        if (transferIndex > max_transfer_index)
        {

            max_transfer_index = transferIndex;
            transferIndexesSize.clear();
            transferIndexesSize.shrink_to_fit();

            traLength = seqs[i].size();
            transferIndexesSize.push_back(i);
        }
        else if (transferIndex == max_transfer_index)
        {

            if (seqs[i].size() < traLength)
            {
                transferIndexesSize.clear();
                transferIndexesSize.shrink_to_fit();
                transferIndexesSize.push_back(i);
                traLength = seqs[i].size();
            }
            if (seqs[i].size() == traLength)
            {

                transferIndexesSize.push_back(i);

            }
        }
    }
    if (transferIndexesSize.size() == 0)
    {
        return -1;
    }
    else
    {
        int rI = transferIndexesSize[rand() % transferIndexesSize.size()];
        return rI;
    }
    return -1;
}
void removePrefixesClassicS(vector<vector<vector<int>>>& CS)//Remove all sequences in CS that are prefixes of other sequences (variant S).
{
    vector<vector<int>> myV;
    for (int asd = 0; asd < CS.size(); asd++)
    {
        for (int awe = 0; awe < CS[asd].size(); awe++)
        {
            myV.push_back(CS[asd][awe]);
        }
    }

    vector<vector<int>> res;

    CS.clear();
    CS.shrink_to_fit();
    CS.push_back(res);
    for (int i = 0; i < myV.size()-1; i++)
    {
        for (int j = i + 1; j < myV.size(); j++)
        {
            if (j == 30)
                int asd = 0;
            if (myV[i].size() >= myV[j].size())
                if (isPrefix(myV[j], myV[i]))
                {
                    myV[j].clear();
                }                 
            if(myV[i].size() < myV[j].size())
                if (isPrefix(myV[i], myV[j]))
                {
                    myV[i].clear();
                    break;
                }
                    
        }
    }
    for (int i = 0; i < myV.size(); i++)
    {
        if (myV[i].size() > 0)
        {
            CS[0].push_back(myV[i]);
        }
    }

}
void removePrefixes(vector<vector<vector<int>>>& CS)//Remove all sequences in CS[0] that are prefixes of other sequences.
{


    {
        if (CS[0].size() > 0)
        {
            for (int k = 0; k < CS[0].size() - 1; k++)
            {
                for (int l = k + 1; l < CS[0].size(); l++)
                {
                    if (CS[0][k].size() > 0 && CS[0][l].size() > 0)
                    {
                        if (CS[0][k].size() < CS[0][l].size() && isPrefix(CS[0][k], CS[0][l]))
                        {
                            CS[0][k].clear();
                            CS[0][k].shrink_to_fit();
                        }
                        else if (CS[0][l].size() < CS[0][k].size() && isPrefix(CS[0][l], CS[0][k]))
                        {
                            CS[0][l].clear();
                            CS[0][l].shrink_to_fit();
                        }
                        else if (CS[0][l].size() == CS[0][k].size() && isPrefix(CS[0][l], CS[0][k]))
                        {
                            CS[0][l].clear();
                            CS[0][l].shrink_to_fit();
                        }
                    }
                }
            }
        }
    }
}
void removePrefixesClassic(vector<vector<vector<int>>>& CS)//Classic version of prefix removal for CS[0].
{


    {

        {
            for (int k = 0; k < CS[0].size() - 1; k++)
            {
                for (int l = k + 1; l < CS[0].size(); l++)
                {
                    if (CS[0][k].size() > 0 && CS[0][l].size() > 0)
                    {
                        if (CS[0][k].size() < CS[0][l].size() && isPrefix(CS[0][k], CS[0][l]))
                        {
                            CS[0][k].clear();
                            CS[0][k].shrink_to_fit();
                        }
                        else if (CS[0][l].size() < CS[0][k].size() && isPrefix(CS[0][l], CS[0][k]))
                        {
                            CS[0][l].clear();
                            CS[0][l].shrink_to_fit();
                        }
                        else if (CS[0][l].size() == CS[0][k].size() && isPrefix(CS[0][l], CS[0][k]))
                        {
                            CS[0][l].clear();
                            CS[0][l].shrink_to_fit();
                        }
                    }
                }
            }
        }
    }
}
void RemoveRepeatations(vector<vector<vector<int>>>& CS)//Remove exact duplicate sequences across CS layers.
{
    for (int i = 0; i < CS.size() - 1; i++)
    {
        for (int j = i + 1; j < CS.size(); j++)
        {
            for (int k = 0; k < CS[i].size(); k++)
            {
                for (int l = 0; l < CS[j].size(); l++)
                {
                    if (CS[i][k].size() > 0 && CS[j][l].size() > 0)
                    {
                        if (CS[j][l].size() == CS[i][k].size() && isPrefix(CS[j][l], CS[i][k]))
                        {
                            CS[j][l].clear();
                            CS[j][l].shrink_to_fit();
                        }
                    }
                }
            }
        }
    }

}
void clearCS(vector<vector<vector<int>>>& CS)//Flatten and clear CS while preserving non-empty sequences.
{
    vector<vector<int>> WSet;
    for (int i = 0; i < CS.size(); i++)
    {
        for (int j = 0; j < CS[i].size(); j++)
        {
            if (CS[i][j].size() > 0)
            {
                WSet.push_back(CS[i][j]);
                CS[i][j].clear();
                CS[i][j].shrink_to_fit();
            }
        }
    }
    CS.clear();
    CS.shrink_to_fit();
    CS.push_back(WSet);

}
void clearCSClassic(vector<vector<vector<int>>>& CS)//Flatten and clear CS[0] while preserving non-empty sequences (classic version).
{
    vector<vector<int>> WSet;

    {
        for (int j = 0; j < CS[0].size(); j++)
        {
            if (CS[0][j].size() > 0)
            {
                WSet.push_back(CS[0][j]);
                CS[0][j].clear();
                CS[0][j].shrink_to_fit();
            }
        }
    }
    CS.clear();
    CS.shrink_to_fit();
    CS.push_back(WSet);

}
void clearCSClassicS(vector<vector<vector<int>>>& CS)//Same as clearCSClassic but applied differently (S variant).
{
    vector<vector<int>> WSet;

    {
        for (int j = 0; j < CS[0].size(); j++)
        {
            if (CS[0][j].size() > 0)
            {
                WSet.push_back(CS[0][j]);
                CS[0][j].clear();
                CS[0][j].shrink_to_fit();
            }
        }
    }
    CS.clear();
    CS.shrink_to_fit();
    CS.push_back(WSet);

}



bool checkPath()//Check whether all transitions in the FSM have been visited.
{
    for (int i = 0; i < edges.size(); i++)
    {
        for (int j = 0; j < edges[i].size(); j++)
        {
            if (edges[i][j].visited == false)
                return false;

        }
    }
    return true;
}
int getFreeSlotCount(int source)//Count the number of unvisited transitions from a given source state.
{
    int fc = 0;
    for (int i = 0; i < edges[0].size(); i++)
    {
        if (edges[source][i].visited == false)
            fc++;
    }
    return fc;
}
int getNextUnvisitedState(int s)//Find the next unvisited state (not equal to current).
{
    for (int i = 0; i < edges.size(); i++)
    {
        for (int j = 0; j < edges[i].size(); j++)
        {
            if (i != s && edges[i][j].visited == false)
                return i;
        }
    }
}
void addIfNotExists(vector<pair<int, vector<int>>>& l, pair<int, vector<int>>& r)//Add pair r to list l if it doesn't already exist.
{
    if (l.size() == 0)
        l.push_back(r);
    else
    {
        for (int i = 0; i < l.size(); i++)
        {
            if (l[i].first == r.first)
                return;
        }
        l.push_back(r);
    }

}
vector<int> recSP(vector<pair<int, vector<int>>>& n, int d)//Recursive helper to find a shortest path from one state to another.
{
    vector<pair<int, vector<int>>> deepen;
    for (int i = 0; i < n.size(); i++)
    {
        for (int j = 0; j < nextStates[0].size(); j++)
        {
            if (nextStates[n[i].first][j] == d)
            {
                n[i].second.push_back(j);
                return n[i].second;
            }
            pair<int, vector<int>> item;
            item.first = nextStates[n[i].first][j];
            item.second = n[i].second;
            item.second.push_back(j);
            if (item.second.size() > FSMStates)
                return vector<int>();
            addIfNotExists(deepen, item);
        }
    }
    if (deepen.size() > nextStates.size())
        return vector<int>();
    else
        return recSP(deepen, d);
}
vector<int> getSP(int s, int d)//Wrapper to get shortest path using BFS-like recursion.
{
    pair<int, vector<int>> p;
    p.first = s;

    vector<pair<int, vector<int>>> t;
    t.push_back(p);
    return recSP(t, d);

}
vector<int> BFS(int s, int d)//Basic BFS to find shortest input sequence from state s to d.
{
    if (s == d)
        return vector<int>();

    vector<int> temp;
    vector<int> sp;
    int currentState = s;
    bool found = false;
    while (!found)
    {

        int input = nextStates[0].size();
        for (int i = 0; i < input; i++)
        {
            if (nextStates[currentState][i] == d)
            {
                temp.push_back(i);
                found = true;
                return temp;
            }
            else
            {
                temp.push_back(i);
                myQ.push(pair<int, vector<int>>(nextStates[currentState][i], temp));
            }
            temp = vector<int>();
            temp = sp;
        }

        pair<int, vector<int>> r = myQ.front();
        myQ.pop();
        vector<int> a;
        sp = r.second;
        temp = sp;
        currentState = r.first;
        if (temp.size() > nextStates.size() - 1)
            return vector<int>();
    }
    return vector<int>();
}
vector<int> getShortestPathFromSource(int s)//Use BFS to find shortest path from source to an unvisited state.
{
    int d = getNextUnvisitedState(s);
    return getSP(s, d);
}

int selectNextInput(int s)//Select the next input for Euler path construction, prioritizing unvisited transitions.
{
    int input = -1;
    int mx_free_slot = INT_MIN;
    vector<int> fss;
    for (int i = 0; i < edges[s].size(); i++)
    {
        if (edges[s][i].visited == false)
        {
            int freeSpace = getFreeSlotCount(edges[s][i].dest);
            if (freeSpace > mx_free_slot)
            {
                fss.clear();
                fss.push_back(i);
                mx_free_slot = freeSpace;
            }
            else if (freeSpace == mx_free_slot)
            {
                mx_free_slot = freeSpace;
                fss.push_back(i);
            }
        }

    }
    if (fss.size() == 0)
    {
        return -1;
    }
    else
    {
        int index = fss[0];
        edges[s][index].visited = true;
        return index;
    }
}


//Construct an Eulerian sequence that visits all edges from a given source.
vector<int> constructEuler(int source, vector<vector<int>>& CS)
{
    bool constructed = false;
    int s = source;
    while (constructed == false)
    {
        indexes.push_back(source);
        int index = selectNextInput(source);
        if (index >= 0)
        {
            for (int k = 0; k < CS[index].size(); k++)
            {
                inputsOrder.push_back(CS[index][k]);
            }

            source = edges[source][index].dest;
        }
        else
        {
            vector<int> transfer = getShortestPathFromSource(source);
            transferLengthO += transfer.size();
            if (transfer.size() == 0)
            {
                return vector<int>();
            }

            for (int k = 0; k < transfer.size(); k++)
            {
                source = nextStates[source][transfer[k]];
                inputsOrder.push_back(k);
            }
        }
        constructed = checkPath();

    }



    return inputsOrder;
}
// DFS 
vector<int> DFS(int source)
{  
    if (true)
        return vector<int>(1, -1);
    for (int i = 0; i < edges[source].size(); i++)
    {
        if (edges[source][i].visited)
        {
            ;
        }
        else
        {
            edges[source][i].visited = true;
            sequence.push(i);
            vector<int> seq = DFS(edges[source][i].dest);
            if (seq.size() > 0)
            {
                seq.push_back(i);
                return seq;
            }
            if (seq.size() == 0)
            {
                edges[source][i].visited = false;
                if (sequence.size() > 0)
                    sequence.pop();

               
            }
        }
    }
    return vector<int>();
}

//Calculate reaching sequences from initial state to all others.
vector<vector<int>> calculateReachingSequneces()
{
    vector<vector<int>> seq;
    int states = nextStates.size();
    seq.resize(states);
    seq[0] = vector<int>();
    for (int i = 1; i < states; i++)
    {
        seq[i] = getSP(0, i);
        transferLengthC += seq[i].size();
        if (seq[i].size() < 1)
            cout << "NO PATH" << endl;
        myQ = queue<pair<int, vector<int>>>();
    }
    return seq;
}
//Generate transition identification sequences using reaching and characterizing sets.
vector<vector<int>> calculateTransitionIdentifiationSequences(vector<vector<int>>& seq, vector<vector<vector<int>>>& CS)
{
    vector<vector<int>>tc;
    int states = nextStates.size();
    int inputs = nextStates[0].size();
    for (int i = 0; i < seq.size(); i++)
    {
        for (int j = 0; j < inputs; j++)
        {
            for (int k = 0; k < CS[0].size(); k++)
            {
                vector<int> temp;
                temp = seq[i];

                temp.push_back(j);
                for (int l = 0; l < CS[0][k].size(); l++)
                {
                    temp.push_back(CS[0][k][l]);
                }
                tc.push_back(temp);
                temp.clear();
                temp.shrink_to_fit();
            }
        }
    }
    return tc;
}
//Generate state identification sequences using reaching and characterizing sets.
vector<vector<int>> calculateStateIdentificationSequences(vector<vector<int>>& seq, vector<vector<vector<int>>>& CS)
{
    vector<vector<int>> result;
    for (int i = 0; i < seq.size(); i++)
    {
        for (int j = 0; j < CS[0].size(); j++)
        {
            vector<int> temp = seq[i];
            for (int k = 0; k < CS[0][j].size(); k++)
            {
                temp.push_back(CS[0][j][k]);
            }
            result.push_back(temp);
            temp = vector<int>();
        }
    }
    return result;
}
//Check if Eulerian path exists(i.e., all transitions are balanced).
bool hasEuler()
{
    int states = nextStates.size();
    vector<int> reached;
    reached.resize(states);
    reached.assign(states, edges[0].size());
    for (int j = 0; j < edges.size(); j++)
    {
        for (int i = 0; i < edges[0].size(); i++)
        {
            if (edges[j][i].dest != j)
                reached[edges[j][i].dest]--;
        }
    }
    int ctr = 0;
    for (int j = 0; j < states; j++)
    {
        if (reached[j] > 0 || reached[j] < 0)
            ctr++;
        if(ctr>1)
            return false;
    }

    return true;
}
//Construct an ordered State Identification Sequence from CS using an Euler path.
vector<int> constructOrderedSIS(vector<vector<vector<int>>>& CS, int FSMInputs)
{
    if (CS[0].size() == 0)
    {
        return vector<int>();
    }
    vector<int> sCover;
    int initialState = 0;
    vector<vector<int>> sequence;

    int states = nextStates.size();
    edges.resize(states);
    vector<edge> temp;
    temp.resize(CS[0].size());
    edges.assign(states, temp);

    for (int i = 0; i < states; i++)
    {
        for (int j = 0; j < CS[0].size(); j++)
        {
            int nextState = reachedState(i, CS[0][j]);
            edges[i][j].dest = nextState;
            edges[i][j].visited = false;
        }
    }
    return constructEuler(0, CS[0]);
}
//Merge SIS and TIS, filter them to retain only characterizing sequences.
vector<vector<int>> calculateCharacterisingSet(vector<vector<int>>& sis, vector<vector<int>>& tis)
{
    vector<vector<int>> ttis = sis;
    for (int i = 0; i < tis.size(); i++)
    {
        ttis.push_back(tis[i]);
    }
    vector<vector<int>> W;

    for (int i = 0; i < ttis.size(); i++)
    {
        vector<int> t = ttis[i];
        bool add = true;
        for (int j = 0; j < ttis.size(); j++)
        {
            if (i != j && t.size() < ttis[j].size())
            {
                if (isPrefix(t, ttis[j]))
                {
                    add = false;
                    j = ttis.size();
                }
            }
        }
        if (add)
            W.push_back(ttis[i]);
    }
    return W;
}
//Variant: takes SIS as a single sequence instead of a set.
vector<vector<int>> calculateCharacterisingSet(vector<int>& sis, vector<vector<int>>& tis)
{
    vector<vector<int>> ttis;
    ttis.push_back(sis);
    for (int i = 0; i < tis.size(); i++)
    {
        ttis.push_back(tis[i]);
    }
    vector<vector<int>> W;

    for (int i = 0; i < ttis.size(); i++)
    {
        vector<int> t = ttis[i];
        bool add = true;
        for (int j = 0; j < ttis.size(); j++)
        {
            if (i != j && t.size() < ttis[j].size())
            {
                if (isPrefix(t, ttis[j]))
                {
                    add = false;
                    j = ttis.size();
                }
            }
        }
        if (add)
            W.push_back(ttis[i]);
    }
    return W;
}
//Get total length of all sequences in a vector.
int getLength(vector<vector<int>>& v)
{
    int result = 0;
    for (int i = 0; i < v.size(); i++)
    {
        result += v[i].size();
    }
    return result;
}
//Check if a sequence v2 is a prefix of any in v1.
bool isFamily(vector<vector<int>>& v1, vector<int>& v2)
{

    for (int i = 0; i < v1.size(); i++)
    {
        if (isLinR(v2, v1[i]))
        {
            true;
        }
    }
    return false;
}
//Count how many times a sequence appears as a prefix in other partitions.
int GetOccurrance(int l, vector<int>& myString, vector<vector<vector<int>>>& SSCondesed)
{
    int observed = 0;
    for (int i = 0; i < SSCondesed.size(); i++)
    {
        if (i != l)
        {
            for (int j = 0; j < SSCondesed[i].size(); j++)
            {
                if (SSCondesed[i][j].size() > 0 && isLinR(SSCondesed[i][j], myString))
                    observed++;
            }
        }
    }
    return observed;
}
// Count how many sequences in SSCondesed[0] have the input myString as a prefix up to length l.
int GetPrefix(int l, vector<int>& myString, vector<vector<vector<int>>>& SSCondesed)
{
    int observed = 0;


    for (int j = 0; j < SSCondesed[0].size(); j++)
    {
        if (isLinR(myString, SSCondesed[0][j]))
            observed++;
    }


    return observed;
}
// Get a random input index that distinguishes the observed outputs of states i and j.
// Returns -1 if no such input exists.
int getSepInfo(int i, int j)
{
    vector<int> sepIn;
    for (int k = 0; k < nextStates[0].size(); k++)
    {
        if (observedOutputs[i][k] != observedOutputs[j][k])
            sepIn.push_back(k);
    }
    if (sepIn.size() == 0)
        return -1;
    else
        return sepIn[rand() % sepIn.size()];
}
// Determine if the given input sequence `seq` distinguishes state i from state j.
// Returns true if the sequence leads to different outputs or transitions.
bool separates(int i, int j, vector<int>& seq)
{
    int o1 = -1;
    int o2 = -1;
    int index = 0;
    while (o1 == o2)
    {
        if (index >= seq.size())
            return false;
        o1 = observedOutputs[i][seq[index]];
        o2 = observedOutputs[j][seq[index]];
        i = nextStates[i][seq[index]];
        j = nextStates[j][seq[index]];

        index++;

    }
    return true;
}
// Check if a single input `in` distinguishes state s1 from state s2 based on observed output.

bool separates(int in, int s1, int s2)
{
    if (observedOutputs[s1][in] != observedOutputs[s2][in])
        return true;
}


int main()
{
    // Enable debug heap allocations and automatic memory leak check at program exit (Visual Studio-specific).
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

    // Variable declarations and initializations.
    int st;
    srand(time(0));  // Seed the random number generator with the current time.

    // FSM (Finite State Machine) configuration parameters.
    FSMid = 0, FSMStates = 0, FSMSize = 0, FSMInputs = 0, FSMOutputs = 0, dv = 0;

    // Temporary variables used for processing FSM transitions and I/O.
    int currentState = 0, nextState = 0, currentInput = 0, currentOutput = 0;
    int counter = 0;        // General-purpose counter
    int it = 0;             // Iterator or loop control variable

    // Prompt user for FSM input file name.
    cout << "Enter FSM file name" << endl;
    string FSM_File_Name;
    cin >> FSM_File_Name;

    // Construct filenames for statistics and sequences output.
    string StatF = FSM_File_Name + "Stat.txt";
    string SeqF = FSM_File_Name + "Seq.txt";

    // Open output files.
    // `writeStatFile` appends to statistics file, preserving previous data.
    // `writeSeqFile` creates a new file for writing sequences.
    ofstream writeStatFile(StatF, ios::app);
    ofstream writeSeqFile(SeqF);

    // Append file extension for FSM description file.
    FSM_File_Name += ".txt";

    // Open the FSM input file for reading.
    ifstream reader(FSM_File_Name);

    // Reset FSM parameters and processing variables again before reading input.
    FSMid = 0, FSMStates = 0, FSMSize = 0, FSMInputs = 0, FSMOutputs = 0, dv = 0;
    currentState = 0, nextState = 0, currentInput = 0, currentOutput = 0;

    // Temporary container and flags for FSM parsing.
    vector<int> temp;       // Used for temporary data storage during FSM parsing.
    int ccc = -1;           // A control flag.
    bool printed = true;    // Flag to manage conditional printing/output.
    while (!reader.eof())
    {
        ccc++;                           // Increment a control counter used for logging progress.
        if (ccc == 100)                 // Reset the counter every 100 FSMs processed.
            ccc = 0;

        int Fi, Fs, Fz, Fn, Fo;         // Temporary holders for FSM metadata.

        temp.clear();                   // Clear temporary vector used for default row initialization.
        temp.shrink_to_fit();           // Force memory deallocation for temp (optimization).

        // Read FSM metadata from input file:
        // Format expected: ID, number of states, number of transitions, inputs, outputs, "dv" (likely a flag or version).
        reader >> FSMid >> FSMStates >> FSMSize >> FSMInputs >> FSMOutputs >> dv;

        // Output current progress to console for monitoring.
        cout << "Processing ctr : " << ccc << " with id " << FSMid << " ";
        printed = false;                // Set printed flag to false (likely used later to control conditional printing).

        // Store metadata in separate variables for later use.
        Fi = FSMid;
        Fs = FSMStates;
        Fz = FSMSize;
        Fn = FSMInputs;
        Fo = FSMOutputs;

        // Initialize `temp` vector to default input size with -1s, used as row prototype.
        temp.resize(FSMInputs);
        temp.assign(FSMInputs, -1);

        // Prepare and initialize the `nextStates` matrix [states x inputs].
        nextStates.clear();
        nextStates.shrink_to_fit();
        nextStates.resize(FSMStates);
        nextStates.assign(FSMStates, temp);

        // Prepare and initialize the `observedOutputs` matrix [states x inputs].
        observedOutputs.clear();
        observedOutputs.shrink_to_fit();
        observedOutputs.resize(FSMStates);
        observedOutputs.assign(FSMStates, temp);

        // Compute the number of unique state pairs (used for distinguishing states).
        pairCount = FSMStates * (FSMStates - 1) / 2;

        // Clean up and reinitialize separating sequences.
        for (int sS = 0; sS < separatingSequences.size(); sS++)
        {
            separatingSequences[sS].clear();          // Clear existing sequence for each pair.
            separatingSequences[sS].shrink_to_fit();  // Release memory.
        }
        separatingSequences.clear();                  // Clear the whole container.
        separatingSequences.shrink_to_fit();          // Free memory.
        separatingSequences.resize(pairCount);        // Resize for new FSM state pairs.

        // Begin reading FSM transitions from the file.
        int counter = 0;
        bool partial = false;         // Flag to track if any transition line is partially malformed.

        while (counter < FSMSize)     // Read FSMSize number of transitions.
        {
            char input;               // Input is read as a character.
            reader >> currentState >> nextState >> input >> currentOutput;

            currentInput = int(input - 97);   // Convert 'a', 'b', ... to 0, 1, ... (ASCII offset).
            currentState--;                   // Convert to 0-based indexing.
            nextState--;

            // If either state is invalid (negative), mark FSM as partially defined.
            if (currentState < 0 || nextState < 0)
            {
                partial = true;
            }
            else
            {
                // Store valid transition and output.
                nextStates[currentState][currentInput] = nextState;
                observedOutputs[currentState][currentInput] = currentOutput;
            }

            counter++;
        }

        // Sanity check: alert if FSM has invalid state count.
        if (FSMStates < 0)
        {
            cerr << "here";    // Basic debug message — could be expanded for clarity.
        }

        else
        {

            // Declare data structures for characterization sets, witnesses, and sequence analysis.
            vector<vector<vector<int>>> CSO;  // Characterizing Set - Ordered
            vector<vector<vector<int>>> CSN;  // Characterizing Set - Non-deterministic or Intermediate
            vector<vector<vector<int>>> CSM;  // Characterizing Set - Minimal

            vector<vector<int>> WN;  // Witnesses for Non-deterministic sequences
            vector<vector<int>> WO;  // Witnesses for Ordered sequences
            vector<vector<int>> WM;  // Witnesses for Minimal sequences

            vector<vector<int>> reachingSequences;  // Input sequences reaching each state
            vector<vector<int>> tis;                // Transition Identification Sequences
            vector<vector<int>> sis;                // State Identification Sequences

            // Initialize active states if FSM is not partial (fully defined).
            vector<int> activeStates;
            for (int l = 0; l < FSMStates && !partial; l++)
                activeStates.push_back(l);

            // Variables for separating sequence classes and processing.
            vector<pair<int, int>> Origin; // Holds pairs of distinguishable states.
            vector<pair<vector<pair<int, int>>, vector<int>>> classes; // Holds state pairs with their separating input sequences.
            vector<vector<vector<int>>> CS; // Temporary storage for separating sequences.
            CS.resize(FSMStates);

            // Start timing the ordered characterization process.
            auto start_T_time = std::chrono::high_resolution_clock::now();
            transferLengthO = 0;

            // ============== ORDERED CHARACTERIZATION PROCESS STARTS ============== //
            vector<int> inputSequence;
            for (int l = 0; l < FSMInputs && !partial; l++)
            {
                // Check which state pairs are separated by this input.
                Origin = separates(activeStates, l);

                inputSequence.clear(); // Reuse vector.
                inputSequence.shrink_to_fit();
                inputSequence.push_back(l);

                // Save separating pairs with input sequence.
                classes.push_back(pair<vector<pair<int, int>>, vector<int>>(Origin, inputSequence));
            }

            inputSequence.clear(); // Clean up after initialization loop.
            inputSequence.shrink_to_fit();

            int itC = 0;

            // Expand separating input sequences while state pairs remain unseparated.
            for (int l = 0; pairCount > 0 && !partial; l++)
            {
                itC++;
                for (int m = 0; m < FSMInputs; m++)
                {
                    // Extend each input sequence in the class by one more input and check again.
                    classes.push_back(incrementAClass(classes[l], m));
                }
            }

            // If state pairs are still unseparated, FSM is non-minimal or ambiguous.
            if (pairCount > 0 && !partial)
            {
                // Clear all separating sequence data since analysis cannot continue reliably.
                separatingSequences.clear();
                separatingSequences.shrink_to_fit();

                Origin.clear();
                CS.clear();
                CS.shrink_to_fit();
                classes.clear();
                classes.shrink_to_fit();
                activeStates.clear();
                activeStates.shrink_to_fit();
            }
            else // FSM is minimal; continue with building characterization sets.
            {
                vector<vector<vector<int>>> SSCondesed; // Condensed separating sequences.
                SSCondesed.resize(1);

                // From separatingSequences, pick most effective sequences for each pair.
                for (int l = 0; l < separatingSequences.size(); l++)
                {
                    int t = getOrderingSequence(separatingSequences[l]); // Get index of ordered sequence.
                    if (t == -1)
                    {
                        break; // Stop if no valid ordered sequence.
                    }
                    else
                    {
                        SSCondesed[0].push_back(separatingSequences[l][t]); // Store the selected one.
                    }
                }

                // Remove redundant prefixes from sequences (classic reduction step).
                removePrefixesClassicS(SSCondesed);

                bool end = false;
                vector<int> oSIS;      // Ordered State Identification Sequence.
                int trialCounter = 0;
                int trialMX = 400;     // Max trials for constructing SIS.
                bool fail = false;

                oSIS = constructOrderedSIS(SSCondesed, FSMInputs); // Generate ordered SIS from condensed SS.

                if (!fail)
                {
                    separatingSequences.clear();
                    separatingSequences.shrink_to_fit();

                    if (oSIS.size() > 0) // If we successfully found a valid SIS.
                    {
                        reachingSequences = calculateReachingSequneces(); // Get input paths to all states.
                        tis = calculateTransitionIdentifiationSequences(reachingSequences, SSCondesed); // TIS: distinguish transitions.

                        // Calculate total size of TIS and SIS for reporting.
                        int OtisS = 0;
                        for (int fg = 0; fg < tis.size(); fg++)
                            OtisS += tis[fg].size();
                        int OsisS = oSIS.size();

                        int otiscount = tis.size(); // Number of transitions identified.
                        int osiscount = 1;          // Count of SIS (just one here).

                        WO = calculateCharacterisingSet(oSIS, tis); // Final characterization set.

                        double memoryO = getUsedMemoryMB(); // Report memory usage.

                        // Measure total ordered processing time.
                        auto end_T_time = std::chrono::high_resolution_clock::now();
                        auto OTime_T2 = end_T_time - start_T_time;
                        float t_T1 = std::chrono::duration_cast<std::chrono::milliseconds>(OTime_T2).count();

                        // ============== ORDERED CHARACTERIZATION ENDS ============== //
                                   
                        // Clear all previous data structures to prepare for a new FSM minimization process.
                        separatingSequences.clear();  // Clear and shrink to free memory.
                        separatingSequences.shrink_to_fit();
                        pairCount = FSMStates * (FSMStates - 1) / 2;  // Initialize pair count for FSM states.
                        separatingSequences.resize(pairCount);  // Resize to hold state pairs.

                        reachingSequences.clear();  // Clear reaching sequences.
                        reachingSequences.shrink_to_fit();

                        tis.clear();  // Clear transition identification sequences.
                        tis.shrink_to_fit();

                        sis.clear();  // Clear state identification sequences.
                        sis.shrink_to_fit();

                        oSIS.clear();  // Clear ordered state identification sequences.
                        oSIS.shrink_to_fit();

                        CS.clear();  // Clear characterization sets.
                        CS.shrink_to_fit();

                        activeStates.clear();  // Clear active states.
                        activeStates.shrink_to_fit();

                        Origin.clear();  // Clear state origin information.
                        Origin.shrink_to_fit();

                        classes.clear();  // Clear state classes.
                        classes.shrink_to_fit();

                        // Initialize active states with all FSM states, unless FSM is partial.
                        for (int l = 0; l < FSMStates && !partial; l++)
                            activeStates.push_back(l);

                        // Resize the CS vector to store separating sequences for each state.
                        CS.resize(FSMStates);

                        // Start timing the classic FSM minimization process.
                        auto start_classic_T_time = std::chrono::high_resolution_clock::now();

                        // ============== CLASSIC STARTS HERE ============== //

                        vector<pair<pair<int, int>, vector<int>>> splittingInfo;  // Stores splitting information for state pairs.
                        for (int l = 0; l < FSMStates - 1; l++)  // Loop through all state pairs (l, m).
                        {
                            for (int m = l + 1; m < FSMStates; m++)
                            {
                                int input = getSepInfo(l, m);  // Get separating input for states l and m.
                                if (input >= 0)  // If there is an input that separates these states.
                                {
                                    pair<int, int> st;
                                    st.first = l;  // First state in pair.
                                    st.second = m; // Second state in pair.

                                    pair<pair<int, int>, vector<int>> sep;
                                    sep.first = st;  // Store state pair.

                                    vector<int> in;  // Vector to store separating input sequence.
                                    in.push_back(input);  // Add the separating input.

                                    CS[0].push_back(in);  // Add input sequence to the characterization set.

                                    sep.second = in;  // Store input sequence in splitting information.
                                    splittingInfo.push_back(sep);

                                    // Find index of the state pair in the separatingSequences array.
                                    int index = returnPairIndex(l, m, observedOutputs.size());
                                    if (separatingSequences[index].size() == 0)
                                    {
                                        pairCount--;  // Decrease pair count as this pair is now separated.
                                        separatingSequences[index].push_back(in);  // Add separating sequence.
                                    }
                                }
                            }
                        }

                        // Keep iterating until all state pairs are separated.
                        while (pairCount > 0)
                        {
                            for (int l = 0; l < FSMStates - 1; l++)  // Loop through state pairs again.
                            {
                                for (int m = l + 1; m < FSMStates; m++)
                                {
                                    int index = returnPairIndex(l, m, observedOutputs.size());  // Get pair index.
                                    if (separatingSequences[index].size() == 0)  // If this pair is still not separated.
                                    {
                                        vector<int> divider;  // Store inputs that can separate the pair.

                                        // Try all inputs to find one that separates states.
                                        for (int n = 0; n < FSMInputs; n++)
                                        {
                                            int ns1 = nextStates[l][n];  // Next state for l on input n.
                                            int ns2 = nextStates[m][n];  // Next state for m on input n.

                                            // Swap if necessary to maintain order.
                                            if (ns1 > ns2)
                                            {
                                                int t = ns1;
                                                ns1 = ns2;
                                                ns2 = t;
                                            }

                                            // If next states are different, look for a separator.
                                            if (ns1 != ns2)
                                            {
                                                int nIndex = returnPairIndex(ns1, ns2, nextStates.size());
                                                if (separatingSequences[nIndex].size() > 0)
                                                {
                                                    divider.push_back(n);  // Add input to separator list.
                                                }
                                            }
                                        }

                                        if (divider.size() > 0)  // If a separating input was found.
                                        {
                                            int input = divider[rand() % divider.size()];  // Randomly select a separating input.
                                            int ns1 = nextStates[l][input];  // Next state for l after input.
                                            int ns2 = nextStates[m][input];  // Next state for m after input.

                                            // Swap if necessary to maintain order.
                                            if (ns1 > ns2)
                                            {
                                                int t = ns1;
                                                ns1 = ns2;
                                                ns2 = t;
                                            }

                                            if (ns1 != ns2)  // If next states are different.
                                            {
                                                int nIndex = returnPairIndex(ns1, ns2, nextStates.size());

                                                vector<int> prevSeq;
                                                prevSeq.push_back(input);  // Start the sequence with the separating input.
                                                for (int o = 0; o < separatingSequences[nIndex][0].size(); o++)
                                                {
                                                    prevSeq.push_back(separatingSequences[nIndex][0][o]);  // Append existing sequence.
                                                }

                                                pairCount--;  // Decrease pair count.
                                                separatingSequences[index].push_back(prevSeq);  // Add to separating sequences.
                                                CS[0].push_back(prevSeq);  // Add to characterization set.

                                                pair<pair<int, int>, vector<int>> sep;
                                                sep.first = { l, m };  // Store state pair.
                                                sep.second = prevSeq;  // Store input sequence for separation.
                                                splittingInfo.push_back(sep);
                                            }
                                        }

                                        divider.clear();  // Clear the divider list.
                                        divider.shrink_to_fit();
                                    }
                                }
                            }
                        }

                        // After all pairs are separated, remove prefixes from sequences.
                        removePrefixesClassic(CS);  // Remove redundant prefixes from characterization set.

                        // Clear and finalize CS after minimization.
                        clearCSClassic(CS);

                        // Calculate reaching sequences (sequences that lead to each state).
                        reachingSequences = calculateReachingSequneces();

                        // Finalize data structures for FSM minimization results.
                        CSN = CS;  // Store final characterization set.
                        tis = calculateTransitionIdentifiationSequences(reachingSequences, CS);  // Calculate transition sequences.
                        sis = calculateStateIdentificationSequences(reachingSequences, CS);  // Calculate state identification sequences.

                        // Calculate statistics on the sequences for reporting.
                        int NtisS = 0;
                        for (int fg = 0; fg < tis.size(); fg++)
                            NtisS += tis[fg].size();  // Total number of transition identification sequences.

                        int NsisS = 0;
                        for (int fg = 0; fg < sis.size(); fg++)
                            NsisS += sis[fg].size();  // Total number of state identification sequences.

                        int ntiscount = tis.size();  // Count of transition identification sequences.
                        int nsiscount = sis.size();  // Count of state identification sequences.

                        // Calculate characterization set (WN) based on SIS and TIS.
                        WN = calculateCharacterisingSet(sis, tis);

                        // Report memory usage during the minimization process.
                        double memoryC = getUsedMemoryMB();

                        // End timing the classic FSM minimization process.
                        auto end_classic_T_time = std::chrono::high_resolution_clock::now();
                        auto classicTime_T2 = end_classic_T_time - start_classic_T_time;
                        float t_T2 = std::chrono::duration_cast<std::chrono::milliseconds>(classicTime_T2).count();

                        // Transfer length variables for future use.
                        transferLengthM = transferLengthC;
                        transferLengthC = 0;

                        // ============== CLASSIC ENDS HERE ============== //

                                        
                       // Clearing and shrinking all previous data structures to free memory and prepare for a new FSM minimization process.
                        separatingSequences.clear();  // Clears the separating sequences.
                        separatingSequences.shrink_to_fit();  // Shrinks the capacity to fit the current size.
                        reachingSequences.clear();  // Clears the reaching sequences.
                        reachingSequences.shrink_to_fit();  // Shrinks the capacity to fit the current size.
                        tis.clear();  // Clears the transition identification sequences.
                        tis.shrink_to_fit();  // Shrinks the capacity to fit the current size.
                        sis.clear();  // Clears the state identification sequences.
                        sis.shrink_to_fit();  // Shrinks the capacity to fit the current size.
                        oSIS.clear();  // Clears the ordered state identification sequences.
                        oSIS.shrink_to_fit();  // Shrinks the capacity to fit the current size.
                        CS.clear();  // Clears the characterization sets.
                        CS.shrink_to_fit();  // Shrinks the capacity to fit the current size.
                        activeStates.clear();  // Clears the list of active states.
                        activeStates.shrink_to_fit();  // Shrinks the capacity to fit the current size.
                        Origin.clear();  // Clears the state origins.
                        Origin.shrink_to_fit();  // Shrinks the capacity to fit the current size.
                        classes.clear();  // Clears the state classes.
                        classes.shrink_to_fit();  // Shrinks the capacity to fit the current size.
                        splittingInfo.clear();  // Clears the splitting information.
                        splittingInfo.shrink_to_fit();  // Shrinks the capacity to fit the current size.

                        // Initialize the pair count for FSM states (for possible state pair separations).
                        pairCount = FSMStates * (FSMStates - 1) / 2;  // Total number of state pairs to be checked.
                        separatingSequences.resize(pairCount);  // Resizes separatingSequences to accommodate the pairs.


                        // Initialize active states with all FSM states (from 0 to FSMStates-1).
                        for (int l = 0; l < FSMStates; l++)
                            activeStates.push_back(l);

                        // Resize the characterization set (CS) to hold FSM states.
                        CS.resize(FSMStates);

                        // Start timing the FSM minimization process to track the execution time.
                        auto start_min_T_time = std::chrono::high_resolution_clock::now();

                        // ============== FSM MWA STARTS HERE ============== //

                        // Iterate over all pairs of FSM states (l, m).
                        for (int l = 0; l < FSMStates - 1; l++)  // Loop over all states l (except the last one).
                        {
                            for (int m = l + 1; m < FSMStates; m++)  // Loop over all states m (from l+1 to FSMStates).
                            {
                                int input = getSepInfo(l, m);  // Get the separating input for the state pair (l, m).
                                if (input >= 0)  // If a separating input exists (input >= 0 means it's a valid separator).
                                {
                                    pair<int, int> st;
                                    st.first = l;  // Store the first state of the pair.
                                    st.second = m;  // Store the second state of the pair.

                                    pair<pair<int, int>, vector<int>> sep;
                                    sep.first = st;  // Store the state pair in splitting information.

                                    vector<int> in;
                                    in.push_back(input);  // Push the separating input into the vector.

                                    CS[0].push_back(in);  // Add the input sequence to the characterization set (CS).

                                    sep.second = in;  // Store the separating input sequence in splittingInfo.
                                    splittingInfo.push_back(sep);  // Add the pair to the splitting information.

                                    // Check if the state pair (l, m) is already separated, and if not, add it.
                                    int index = returnPairIndex(l, m, observedOutputs.size());
                                    if (separatingSequences[index].size() == 0)
                                    {
                                        pairCount--;  // Decrease the number of pairs to be separated.
                                        separatingSequences[index].push_back(in);  // Add the separating sequence for this pair.
                                    }

                                    // Check for all other state pairs (n, o) and separate them if the separating input divides them.
                                    for (int n = l; n < FSMStates - 1; n++)
                                    {
                                        for (int o = m + 1; o < FSMStates; o++)
                                        {
                                            if (separates(input, n, o))  // If the input separates states n and o.
                                            {
                                                pair<int, int> sts;
                                                sts.first = n;  // Store state n.
                                                sts.second = o;  // Store state o.

                                                pair<pair<int, int>, vector<int>> seps;
                                                seps.first = st;  // Store the state pair in splitting information.
                                                seps.second = in;  // Store the separating input sequence.

                                                splittingInfo.push_back(seps);  // Add the separation info to the splitting list.

                                                int indexs = returnPairIndex(n, o, observedOutputs.size());
                                                if (separatingSequences[indexs].size() == 0)
                                                {
                                                    pairCount--;  // Decrease the pair count.
                                                    separatingSequences[indexs].push_back(in);  // Add the separating sequence for this pair.
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        // Iterate until all pairs have been separated.
                        while (pairCount > 0)
                        {
                            for (int l = 0; l < FSMStates - 1; l++)  // Loop over all states l.
                            {
                                for (int m = l + 1; m < FSMStates; m++)  // Loop over all states m.
                                {
                                    int index = returnPairIndex(l, m, observedOutputs.size());  // Get the index for the pair (l, m).
                                    if (separatingSequences[index].size() == 0)  // If no separating sequence has been found for this pair.
                                    {
                                        vector<int> divider;  // Store possible separating inputs.

                                        // Loop over all FSM inputs to find which input separates the pair.
                                        for (int n = 0; n < FSMInputs; n++)
                                        {
                                            int ns1 = nextStates[l][n];  // Next state for l on input n.
                                            int ns2 = nextStates[m][n];  // Next state for m on input n.

                                            // Ensure states are ordered (swap if necessary).
                                            if (ns1 > ns2)
                                            {
                                                int t = ns1;
                                                ns1 = ns2;
                                                ns2 = t;
                                            }

                                            // If the next states differ, check for a separating input.
                                            if (ns1 != ns2)
                                            {
                                                int nIndex = returnPairIndex(ns1, ns2, nextStates.size());
                                                if (separatingSequences[nIndex].size() > 0)
                                                {
                                                    divider.push_back(n);  // Add the input to the divider list.
                                                }
                                            }
                                        }

                                        if (divider.size() > 0)  // If a separating input is found.
                                        {
                                            int input = divider[rand() % divider.size()];  // Randomly select a separating input.
                                            int ns1 = nextStates[l][input];  // Next state for l after input.
                                            int ns2 = nextStates[m][input];  // Next state for m after input.

                                            // Swap if necessary to maintain order.
                                            if (ns1 > ns2)
                                            {
                                                int t = ns1;
                                                ns1 = ns2;
                                                ns2 = t;
                                            }

                                            // If next states differ, separate the pair and add the sequence.
                                            if (ns1 != ns2)
                                            {
                                                int nIndex = returnPairIndex(ns1, ns2, nextStates.size());

                                                vector<int> prevSeq;
                                                prevSeq.push_back(input);  // Start the sequence with the separating input.
                                                for (int o = 0; o < separatingSequences[nIndex][0].size(); o++)
                                                {
                                                    prevSeq.push_back(separatingSequences[nIndex][0][o]);  // Append existing sequences.
                                                }

                                                pairCount--;  // Decrease the number of remaining state pairs.
                                                separatingSequences[index].push_back(prevSeq);  // Store the new sequence for this pair.
                                                CS[0].push_back(prevSeq);  // Add the new sequence to the characterization set.

                                                pair<int, int> st;
                                                st.first = l;
                                                st.second = m;

                                                pair< pair<int, int>, vector<int>>sep;
                                                sep.first = st;
                                                sep.second = prevSeq;
                                                splittingInfo.push_back(sep);  // Add the separation information for this pair.
                                            }
                                        }
                                        divider.clear();  // Clear the divider list.
                                        divider.shrink_to_fit();
                                    }
                                }
                            }
                        }

                        // After all pairs are separated, remove redundant prefixes from the characterization set.
                        removePrefixesClassic(CS);  // Remove any prefixes from the sequences.

                        // Clear the characterization set and finalize the FSM minimization process.
                        clearCSClassic(CS);

                        // Calculate the reaching sequences (sequences that lead to each state).
                        reachingSequences = calculateReachingSequneces();

                        // Store the final characterization set (CS) and compute the transition and state identification sequences.
                        CSM = CS;
                        tis = calculateTransitionIdentifiationSequences(reachingSequences, CS);
                        sis = calculateStateIdentificationSequences(reachingSequences, CS);

                        // Calculate the total number of transition and state identification sequences.
                        int MtisS = 0;
                        for (int fg = 0; fg < tis.size(); fg++)
                            MtisS += tis[fg].size();

                        int MsisS = 0;
                        for (int fg = 0; fg < sis.size(); fg++)
                            MsisS += sis[fg].size();

                        // Count the number of transition and state identification sequences.
                        int Mtiscount = tis.size();
                        int Msiscount = sis.size();

                        // Calculate the characterization set using SIS and TIS.
                        WM = calculateCharacterisingSet(sis, tis);

                        // Track memory usage during the minimization process.
                        double memoryM = getUsedMemoryMB();

                        // End timing the FSM minimization process.
                        auto end_min_T_time = std::chrono::high_resolution_clock::now();
                        auto classicTime_T3 = end_min_T_time - start_min_T_time;
                        float t_T3 = std::chrono::duration_cast<std::chrono::milliseconds>(classicTime_T3).count();

                        // ============== FSM MWA ENDS HERE ============== //

                        // Clear and shrink memory for various data structures
                        separatingSequences.clear();  // Clears the separatingSequences vector
                        separatingSequences.shrink_to_fit();  // Shrinks the capacity of separatingSequences to fit its size

                        reachingSequences.clear();  // Clears the reachingSequences vector
                        reachingSequences.shrink_to_fit();  // Shrinks the capacity of reachingSequences to fit its size

                        tis.clear();  // Clears the tis vector
                        tis.shrink_to_fit();  // Shrinks the capacity of tis to fit its size

                        sis.clear();  // Clears the sis vector
                        sis.shrink_to_fit();  // Shrinks the capacity of sis to fit its size

                        oSIS.clear();  // Clears the oSIS vector
                        oSIS.shrink_to_fit();  // Shrinks the capacity of oSIS to fit its size

                        CS.clear();  // Clears the CS vector
                        CS.shrink_to_fit();  // Shrinks the capacity of CS to fit its size

                        activeStates.clear();  // Clears the activeStates vector
                        activeStates.shrink_to_fit();  // Shrinks the capacity of activeStates to fit its size

                        Origin.clear();  // Clears the Origin vector
                        Origin.shrink_to_fit();  // Shrinks the capacity of Origin to fit its size

                        classes.clear();  // Clears the classes vector
                        classes.shrink_to_fit();  // Shrinks the capacity of classes to fit its size

                        splittingInfo.clear();  // Clears the splittingInfo vector
                        splittingInfo.shrink_to_fit();  // Shrinks the capacity of splittingInfo to fit its size

                        // Calculate the lengths of various vectors
                        int WOlength = getLength(WO);  // Get the length of WO
                        int WNlength = getLength(WN);  // Get the length of WN

                        int CSNlength = getLength(CSN[0]);  // Get the length of the first element in CSN
                        int WMlength = getLength(WM);  // Get the length of WM
                        int CSMlength = getLength(CSM[0]);  // Get the length of the first element in CSM

                                        
                        {
                            // Ensure that memory values are positive to avoid errors
                            if (memoryO <= 0) memoryO = 0.0001;  // Adjust memoryO if it's less than or equal to zero
                            if (memoryC <= 0) memoryC = 0.0001;  // Adjust memoryC if it's less than or equal to zero
                            if (memoryM <= 0) memoryM = 0.0001;  // Adjust memoryM if it's less than or equal to zero

                            // Check if certain conditions are met for FSM processing
                            if (ccc == FSMid && FSMid == Fi && Fo == FSMOutputs && Fz == FSMSize && Fn == FSMInputs)
                            {
                                cout << FSMid << " " << FSMStates << " " << FSMSize << " " << FSMInputs << " " << FSMOutputs << " is processed.." << endl;
                                printed = true;  // Indicate that FSM was processed
                            }
                            else
                            {
                                // If FSM is not processed, print a message and write to stat file
                                cout << FSMid << " " << FSMStates << " " << FSMSize << " " << FSMInputs << " " << FSMOutputs << " is not processed.." << endl;
                                writeStatFile << FSMid << " " << FSMStates << " " << FSMInputs << " " << FSMOutputs << " " << FSMSize << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " OWA" << endl;
                                writeStatFile << FSMid << " " << FSMStates << " " << FSMInputs << " " << FSMOutputs << " " << FSMSize << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " CWA" << endl;
                                writeStatFile << FSMid << " " << FSMStates << " " << FSMInputs << " " << FSMOutputs << " " << FSMSize << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " " << -1 << " MWA" << endl;
                            }

                            // If FSM was not processed, log an error message
                            if (!printed)
                                cerr << "here";  // Error message indicating FSM wasn't processed

                            // Write FSM processing results to the stat file
                            writeStatFile << FSMid << " " << FSMStates << " " << FSMInputs << " " << FSMOutputs << " " << FSMSize << " " << t_T1 << " " << memoryO << " " << transferLengthO << " " << 1 << " " << OsisS << " " << otiscount + osiscount << " " << OsisS + OtisS << " OWA" << endl;
                            writeStatFile << FSMid << " " << FSMStates << " " << FSMInputs << " " << FSMOutputs << " " << FSMSize << " " << t_T2 << " " << memoryC << " " << transferLengthM << " " << nsiscount << " " << NsisS << " " << ntiscount + nsiscount << " " << NtisS + NsisS << " CWA" << endl;
                            writeStatFile << FSMid << " " << FSMStates << " " << FSMInputs << " " << FSMOutputs << " " << FSMSize << " " << t_T3 << " " << memoryM << " " << transferLengthC << " " << Msiscount << " " << MsisS << " " << Mtiscount + Msiscount << " " << MtisS + MsisS << " MWA" << endl;


                        }
                                        
                                        
                        // Clear the WO and WN vectors and shrink memory
                        WO.clear();  // Clear the WO vector
                        WO.shrink_to_fit();  // Shrink the capacity of WO

                        WN.clear();  // Clear the WN vector
                        WN.shrink_to_fit();  // Shrink the capacity of WN
                    }
                    // If there's an error, print "verr"
                    else
                    cout << "verr";  // Error message

                    // Clear and shrink CSO and CSN vectors
                    CSO.clear();  // Clear the CSO vector
                    CSN.clear();  // Clear the CSN vector
                    CSO.shrink_to_fit();  // Shrink the capacity of CSO
                    CSN.shrink_to_fit();  // Shrink the capacity of CSN

                }
                                
                // Clear and shrink various other data structures
                CS.clear();  // Clear the CS vector
                CS.shrink_to_fit();  // Shrink the capacity of CS

                activeStates.clear();  // Clear the activeStates vector
                activeStates.shrink_to_fit();  // Shrink the capacity of activeStates

                indexes.clear();  // Clear the indexes vector
                indexes.shrink_to_fit();  // Shrink the capacity of indexes

                edges.clear();  // Clear the edges vector
                edges.shrink_to_fit();  // Shrink the capacity of edges

                oSIS.clear();  // Clear the oSIS vector
                oSIS.shrink_to_fit();  // Shrink the capacity of oSIS

                // Reset the sequence and queue data structures
                mySequence = stack<int>();  // Create a new empty stack for mySequence
                sequence = stack<int>();  // Create a new empty stack for sequence
                myQ = queue<pair<int, vector<int>>>();  // Create a new empty queue for myQ
                sms = vector < vector<bool>>();  // Create a new empty vector for sms

                // Re-initialize other vectors and data structures
                indexes = deque<int>();  // Create a new empty deque for indexes
                inputsOrder = vector<int>();  // Create a new empty vector for inputsOrder

                edges = vector<vector<edge>>();  // Create a new empty vector of edge vectors for edges
                myQ = queue<pair<int, vector<int>>>();  // Create a new empty queue for myQ
                sms = vector<vector<bool>>();  // Create a new empty vector of boolean vectors for sms

                // If FSM was not processed, log an error message
                if (!printed)
                    cerr << "here";  // Error message


            }
        }
                        
    }   
    // Close the files used for writing sequences and statistics
    writeSeqFile.close();  // Close the sequence file
    writeStatFile.close();  // Close the statistics file 
}
