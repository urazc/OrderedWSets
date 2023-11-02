#include<fstream>
#include<iostream>
#include<vector>
#include <stack>
#include <queue>
#include <chrono>
#include <Windows.h>
#include <Psapi.h>
#include <deque>
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


double getUsedMemoryMB()
{
    PROCESS_MEMORY_COUNTERS_EX pmc;
    DWORD ret = GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
    if (ret == 0)
    {
        printf("GetProcessMemoryInfo failed with code %d\n", GetLastError());
    }
    size_t temp2 = pmc.PrivateUsage;
    double mem = static_cast<float>(temp2);
    mem = mem / 1048576.0;

    return mem;
}
unsigned int returnPairIndex(int i, int j, int n)
{
    int result = ((i * n) - (i * (i + 1) / 2)) + (j - i - 1);
    return result;
}
vector<pair<int, int>> separates(vector<int>& s, int x)
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
vector<pair<int, int>> separatesClassic(vector<int>& s, int x)
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
pair< vector<pair<int, int>>, vector<int>> incrementAClass(pair< vector<pair<int, int>>, vector<int>>& myClass, int input)
{
    vector<pair<int, int>> RESULT;
    vector<int> seq;
    seq.push_back(input);

    if (myClass.second.size() == 0)
    {
        return  pair< vector<pair<int, int>>, vector<int>>(vector<pair<int, int>>(), vector<int>());
    }
    for (int i = 0; i < myClass.second.size(); i++)
    {
        seq.push_back(myClass.second[i]);
    }

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
pair< vector<pair<int, int>>, vector<int>> incrementAClassClassic(pair< vector<pair<int, int>>, vector<int>>& myClass, int input)
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
bool isAdded(vector<int>& l, int r)
{
    bool nt = false;
    for (int i = 0; i < l.size(); i++)
    {
        if (l[i] == r)
            return true;

    }
    return false;
}
bool isPrefix(vector<int>& L, vector<int>& R)
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
bool areEqual(vector<int>& L, vector<int>& R)
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
bool isLinR(vector<int>& L, vector<int>& R)
{
    if (L.size() <= R.size())
    {
        if (isPrefix(L, R))
            return true;
    }

    return false;
}
int isIn(vector<int>& L, vector<int>& R)
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
bool hasAsset(vector<vector<int>>& cs, vector<int>& seqs)
{
    for (int i = 0; i < cs.size(); i++)
    {
        int r = isIn(seqs, cs[i]);
        if (r == 1)
            return true;
    }
    return false;
}
int reachedState(int state, vector<int> seq)
{
    for (int i = 0; i < seq.size(); i++)
    {
        state = nextStates[state][seq[i]];
    }
    return state;
}
int getOrderingSequence(vector<vector<int>>& seqs)
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

void removePrefixes(vector<vector<vector<int>>>& CS)
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
void removePrefixesClassic(vector<vector<vector<int>>>& CS)
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
void RemoveRepeatations(vector<vector<vector<int>>>& CS)
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
void clearCS(vector<vector<vector<int>>>& CS)
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
void clearCSClassic(vector<vector<vector<int>>>& CS)
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



bool checkPath()
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
int getFreeSlotCount(int source)
{
    int fc = 0;
    for (int i = 0; i < edges[0].size(); i++)
    {
        if (edges[source][i].visited == false)
            fc++;
    }
    return fc;
}
int getNextUnvisitedState(int s)
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
void addIfNotExists(vector<pair<int, vector<int>>>& l, pair<int, vector<int>>& r)
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
vector<int> recSP(vector<pair<int, vector<int>>>& n, int d)
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
vector<int> getSP(int s, int d)
{
    pair<int, vector<int>> p;
    p.first = s;

    vector<pair<int, vector<int>>> t;
    t.push_back(p);
    return recSP(t, d);

}
vector<int> BFS(int s, int d)
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
vector<int> getShortestPathFromSource(int s)
{
    int d = getNextUnvisitedState(s);
    return getSP(s, d);
}

int selectNextInput(int s)
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

int getLength(vector<vector<int>>& v)
{
    int result = 0;
    for (int i = 0; i < v.size(); i++)
    {
        result += v[i].size();
    }
    return result;
}
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

bool separates(int in, int s1, int s2)
{
    if (observedOutputs[s1][in] != observedOutputs[s2][in])
        return true;
}


int main()
{


    int st;

    
    srand(time(0));
    
        FSMid = 0, FSMStates = 0, FSMSize = 0, FSMInputs = 0, FSMOutputs = 0, dv = 0;
        int currentState = 0, nextState = 0, currentInput = 0, currentOutput = 0;
        
       
        int counter = 0;
        
            int it = 0;
            
         
                    cout << "Enter FSM file name" << endl;
                    string FSM_File_Name;
                    cin >> FSM_File_Name;
                    string StatF = FSM_File_Name + "Stat.txt";
                    string SeqF = FSM_File_Name + "Seq.txt";
                    ofstream writeStatFile(StatF);
                    ofstream writeSeqFile(SeqF);
                    FSM_File_Name+=".txt";
     
                    writeStatFile << "FSM_ID FSMStates FSMInputs FSMOutputs FSMSize W_Time W_Memory SIS_Tr SIS_Size SIS_Cost TS_Size TS_Cost Algorithm" << endl;
                    writeSeqFile << "FSM_ID FSMStates FSMInputs FSMOutputs FSMSize Ordered_WSequences --- Classic_WSequences --- Ordered_CSSequences --- Classic_CSSequences" << endl;

                    ifstream reader(FSM_File_Name);
                    FSMid = 0, FSMStates = 0, FSMSize = 0, FSMInputs = 0, FSMOutputs = 0, dv = 0;
                    currentState = 0, nextState = 0, currentInput = 0, currentOutput = 0;
                    while (!reader.eof())
                    {

        
                        reader >> FSMid >> FSMStates >> FSMSize >> FSMInputs >> FSMOutputs >> dv;

                        vector<int> temp;
                        temp.resize(FSMInputs);
                        temp.assign(FSMInputs, -1);
        
                        nextStates.resize(FSMStates);
                        nextStates.assign(FSMStates, temp); 
        
                        observedOutputs.resize(FSMStates);
                        observedOutputs.assign(FSMStates, temp);
        
                        pairCount = FSMStates * (FSMStates - 1) / 2;
                        separatingSequences.resize(pairCount);

                        int counter = 0;
                        bool partial = false;
                        while (counter < FSMSize)
                        {
                            char input;
                            reader >> currentState >> nextState >> input >> currentOutput;
                           
                            currentInput = int(input - 97);
                            currentState--; nextState--;
                            if (currentState < 0 || nextState<0)
                            {
                                
                                partial = true;
                             
                            }
                            else {
                                nextStates[currentState][currentInput] = nextState;
                                observedOutputs[currentState][currentInput] = currentOutput;
                            }

                            counter++;
                        }
                            
                            vector<int> activeStates;
                            for (int l = 0; l < FSMStates && !partial; l++)
                                activeStates.push_back(l);

                            vector<pair<int, int>> Origin;
                            vector<pair< vector<pair<int, int>>, vector<int>>>classes;
                            vector<vector<vector<int>>> CS;
                            CS.resize(FSMStates);
                            auto start_T_time = std::chrono::high_resolution_clock::now();
                            ///////////////////////ORDERED STARTS HERE/////////////////////////
                            for (int l = 0; l < FSMInputs && !partial; l++)
                            {
                                Origin = separates(activeStates, l);//returns set of states that 
                                vector<int>inputSequence;
                                inputSequence.push_back(l);
                                classes.push_back(pair< vector<pair<int, int>>, vector<int>>(Origin, inputSequence));
                            }
                            int itC = 0;

                            for (int l = 0; itC<5 && pairCount >0 && !partial; l++)
                            {
                                itC++;
                                for (int m = 0; m < FSMInputs; m++)
                                {
                                    classes.push_back(incrementAClass(classes[l], m));
                                }

                            }
                            if (pairCount > 0 && !partial)//if FSM is not minimal
                            {
                                separatingSequences.clear();
                                separatingSequences.shrink_to_fit();

                                Origin.clear();
                                Origin.clear();
                                CS.clear();
                                CS.shrink_to_fit();
                                classes.clear();
                                classes.shrink_to_fit();
                                activeStates.clear();
                                activeStates.shrink_to_fit();
                            }
                            else//otherwise
                            {
                                vector<vector<vector<int>>> SSCondesed;
                                for (int l = 0; l < FSMStates - 1; l++)
                                {
                                    int startingIndex = returnPairIndex(l, l + 1, FSMStates);
                                    int endingIndex = returnPairIndex(l, FSMStates - 1, FSMStates);
                                    vector<vector<vector<int>>> tempS;
                                    tempS.resize(1);
                                    for (int m = startingIndex; m <= endingIndex; m++)
                                    {
                                        int t = getOrderingSequence(separatingSequences[m]);
                                        if (t == -1)
                                        {
                                            m = endingIndex;
                                            l = FSMStates;
                                            break;
                                        }
                                        else
                                        {
                                            tempS[0].push_back(separatingSequences[m][t]);
                                        }

                                    }
                                    removePrefixesClassic(tempS);
                                    clearCSClassic(tempS);
                                    SSCondesed.push_back(tempS[0]);
                                    tempS.clear();
                                    tempS.shrink_to_fit();
                                }
                                bool end = false;
                                vector<int> oSIS;
                                int trialCounter = 0;
                                int trialMX = 40;
                                bool fail = false;
                                while (!end && SSCondesed.size()>0)
                                {

                                    for (int l = 0; l < SSCondesed.size(); l++)
                                    {
                                        if (l == 2)
                                            int asd = 0;
                                        if (SSCondesed[l].size() == 1)
                                        {
                                            CS[0].push_back(SSCondesed[l][0]);
                                        }
                                        else
                                        {

                                            int occurrenceIndex = -1;
                                            int mOccurrence = -1;
                                            int prefixIndex = -1;
                                            int mPrefix = -1;

                                            for (int m = 0; m < SSCondesed[l].size(); m++)
                                            {

                                                vector<int> myString = SSCondesed[l][m];
                                                if (myString.size() > 0)
                                                {
                                                    int occurrance = GetOccurrance(l, myString, SSCondesed);
                                                    int beingPrefix = GetPrefix(l, myString, CS);
                                                    if (occurrance + beingPrefix > mOccurrence)
                                                    {
                                                        mOccurrence = occurrance + beingPrefix;
                                                        occurrenceIndex = m;
                                                    }

                                                }

                                            }
                                            if (occurrenceIndex >= 0)
                                            {
                                                CS[0].push_back(SSCondesed[l][occurrenceIndex]);
                                                SSCondesed[l][occurrenceIndex].clear();
                                                SSCondesed[l][occurrenceIndex].shrink_to_fit();
                                            }

                                            else
                                            {
                                                end = true;
                                            }
                                        }
                                    }
                                    if (CS.size() == 0 || CS[0].size() == 0)
                                        ;
                                    else
                                    {
                                        removePrefixes(CS);
                                        clearCS(CS);
                                        
                                        bool add = true;
                                        for (int ll = 0; ll < FSMStates - 1; ll++)
                                        {
                                            for (int mm = ll + 1; mm < FSMStates; mm++)
                                            {
                                                bool seps = false;
                                                for (int jj = 0; jj < CS[0].size() && !seps; jj++)
                                                {
                                                    seps = separates(ll, mm, CS[0][jj]);
                                                }
                                                if (seps == false)
                                                {
                                                    add = false;
                                                }
                                            }
                                        }
                                        if (add)
                                            oSIS = constructOrderedSIS(CS, FSMInputs);
                                    }


                                    if (oSIS.size() > 0)
                                    {

                                        end = true;
                                    }
                                    else
                                    {
                                        oSIS.clear();
                                        oSIS.shrink_to_fit();
                                        if (CS.size() > 0)
                                        {
                                            CS[0].clear();
                                            CS[0].shrink_to_fit();
                                        }

                                    }
                                    trialCounter++;
                                    if (trialCounter > trialMX)
                                        fail = true;
                                }
                                if (!fail)
                                {
                                
                                    vector<vector<vector<int>>> CSO;
                                    vector<vector<vector<int>>> CSN;
                                    vector<vector<vector<int>>> CSM;


                                    separatingSequences.clear();
                                    separatingSequences.shrink_to_fit();
                                    
                                    if (oSIS.size() > 0)
                                    {
                                        CSO = CS;
                                        vector<vector<int>> WN;
                                        vector<vector<int>> WO;
                                        vector<vector<int>> WM;
                                        vector<vector<int>> reachingSequences = calculateReachingSequneces();
                                        vector<vector<int>> tis = calculateTransitionIdentifiationSequences(reachingSequences, CS);
                                        vector<vector<int>> sis;
                                        
                                        int OtisS = 0;
                                        for (int fg = 0; fg < tis.size(); fg++)
                                            OtisS += tis[fg].size();
                                        int OsisS = oSIS.size();

                                        int otiscount = tis.size();
                                        int osiscount = 1;
                                        WO = calculateCharacterisingSet(oSIS, tis);
                                        double memoryO = getUsedMemoryMB();
                                        auto end_T_time = std::chrono::high_resolution_clock::now();
                                        auto OTime_T2 = end_T_time - start_T_time;
                                        float t_T1 = OTime_T2 / std::chrono::milliseconds(1);
                                        ///////////////////////ORDERED ENDS HERE/////////////////////////

                                        separatingSequences.clear();
                                        separatingSequences.shrink_to_fit();
                                        pairCount = FSMStates * (FSMStates - 1) / 2;
                                        separatingSequences.resize(pairCount);
                                        reachingSequences.clear();
                                        reachingSequences.shrink_to_fit();
                                        tis.clear();
                                        tis.shrink_to_fit();
                                        sis.clear();
                                        sis.shrink_to_fit();
                                        oSIS.clear();
                                        oSIS.shrink_to_fit();


                                        CS.clear();
                                        CS.shrink_to_fit();

                                        activeStates.clear();
                                        activeStates.shrink_to_fit();

                                        for (int l = 0; l < FSMStates; l++)
                                            activeStates.push_back(l);
                                        Origin.clear();
                                        Origin.shrink_to_fit();
                                        classes.clear();
                                        classes.shrink_to_fit();
                                        CS.resize(FSMStates);

                                        auto start_classic_T_time = std::chrono::high_resolution_clock::now();
                                        ///////////////////////CLASSIC STARTS HERE/////////////////////////
                                        vector< pair<pair<int, int>, vector<int>>> splittingInfo;
                                        for (int l = 0; l < FSMStates - 1; l++)
                                        {
                                            for (int m = l + 1; m < FSMStates; m++)
                                            {
                                                int input = getSepInfo(l, m);
                                                if (input >= 0)
                                                {
                                                    pair<int, int> st;
                                                    st.first = l;
                                                    st.second = m;
                                                    pair< pair<int, int>, vector<int>> sep;
                                                    sep.first = st;
                                                    vector<int> in;
                                                    in.push_back(input);
                                                    CS[0].push_back(in);
                                                    sep.second = in;
                                                    splittingInfo.push_back(sep);
                                                    int index = returnPairIndex(l, m, observedOutputs.size());
                                                    if (separatingSequences[index].size() == 0)
                                                    {
                                                        pairCount--;
                                                        separatingSequences[index].push_back(in);
                                                    }

                                                }
                                            }
                                        }
                                        while (pairCount > 0)
                                        {
                                            for (int l = 0; l < FSMStates - 1; l++)
                                            {
                                                for (int m = l + 1; m < FSMStates; m++)
                                                {
                                                    int index = returnPairIndex(l, m, observedOutputs.size());
                                                    if (separatingSequences[index].size() == 0)
                                                    {
                                                        vector<int> divider;
                                                        for (int n = 0; n < FSMInputs; n++)
                                                        {
                                                            int ns1 = nextStates[l][n];
                                                            int ns2 = nextStates[m][n];
                                                            if (ns1 > ns2)
                                                            {
                                                                int t = ns1;
                                                                ns1 = ns2;
                                                                ns2 = t;
                                                            }

                                                            if (ns1 == ns2)
                                                            {
                                                                ;
                                                            }
                                                            else
                                                            {
                                                                int nIndex = returnPairIndex(ns1, ns2, nextStates.size());
                                                                if (separatingSequences[nIndex].size() > 0)
                                                                {
                                                                    divider.push_back(n);
                                                                }
                                                            }

                                                        }
                                                        if (divider.size() > 0)
                                                        {
                                                            int input = divider[rand() % divider.size()];
                                                            int ns1 = nextStates[l][input];
                                                            int ns2 = nextStates[m][input];
                                                            if (ns1 > ns2)
                                                            {
                                                                int t = ns1;
                                                                ns1 = ns2;
                                                                ns2 = t;
                                                            }
                                                            if (ns1 == ns2)
                                                            {
                                                                ;
                                                            }
                                                            else
                                                            {
                                                                int nIndex = returnPairIndex(ns1, ns2, nextStates.size());
                                                                vector<int>prevSeq;
                                                                prevSeq.push_back(input);
                                                                for (int o = 0; o < separatingSequences[nIndex][0].size(); o++)
                                                                {
                                                                    prevSeq.push_back(separatingSequences[nIndex][0][o]);
                                                                }
                                                                if (prevSeq.size() == 0)
                                                                    int ne = 0;
                                                                pairCount--;
                                                                separatingSequences[index].push_back(prevSeq);
                                                                CS[0].push_back(prevSeq);
                                                                pair<int, int> st;
                                                                st.first = l;
                                                                st.second = m;
                                                                pair< pair<int, int>, vector<int>> sep;
                                                                sep.first = st;
                                                                sep.second = prevSeq;
                                                                splittingInfo.push_back(sep);
                                                            }
                                                        }
                                                        else
                                                            ;

                                                        divider.clear();
                                                        divider.shrink_to_fit();
                                                    }
                                                }
                                            }
                                        }

                                        separatingSequences.clear();
                                        splittingInfo.clear();
                                        splittingInfo.shrink_to_fit();
                                        separatingSequences.shrink_to_fit();
                                        removePrefixesClassic(CS);
                                        clearCSClassic(CS);
                                        reachingSequences = calculateReachingSequneces();
                                       
                                        CSN = CS;
                                        tis = calculateTransitionIdentifiationSequences(reachingSequences, CS);
                                        sis = calculateStateIdentificationSequences(reachingSequences, CS);
                                        int NtisS = 0;
                                        for (int fg = 0; fg < tis.size(); fg++)
                                            NtisS += tis[fg].size();
                                        int NsisS = 0;
                                        for (int fg = 0; fg < sis.size(); fg++)
                                            NsisS += sis[fg].size();
                                        int ntiscount = tis.size();
                                        int nsiscount = sis.size();
                                        classes.clear();
                                        classes.shrink_to_fit();
                                        Origin.clear();
                                        Origin.shrink_to_fit();
                                        WN = calculateCharacterisingSet(sis, tis);
                                        double memoryC = getUsedMemoryMB();
                                        auto end_classic_T_time = std::chrono::high_resolution_clock::now();
                                        auto classicTime_T2 = end_classic_T_time - start_classic_T_time;
                                        float t_T2 = classicTime_T2 / std::chrono::milliseconds(1);
                                        ///////////////////////CLASSIC ENDs HERE/////////////////////////

                                        separatingSequences.clear();
                                        separatingSequences.shrink_to_fit();
                                        pairCount = FSMStates * (FSMStates - 1) / 2;
                                        separatingSequences.resize(pairCount);
                                        reachingSequences.clear();
                                        reachingSequences.shrink_to_fit();
                                        tis.clear();
                                        tis.shrink_to_fit();
                                        sis.clear();
                                        sis.shrink_to_fit();
                                        oSIS.clear();
                                        oSIS.shrink_to_fit();


                                        CS.clear();
                                        CS.shrink_to_fit();

                                        activeStates.clear();
                                        activeStates.shrink_to_fit();

                                        for (int l = 0; l < FSMStates; l++)
                                            activeStates.push_back(l);
                                        Origin.clear();
                                        Origin.shrink_to_fit();
                                        classes.clear();
                                        classes.shrink_to_fit();
                                        CS.resize(FSMStates);
                                        splittingInfo.clear();
                                        splittingInfo.shrink_to_fit();
                                        auto start_min_T_time = std::chrono::high_resolution_clock::now();
                                        ///////////////////////Min STARTS HERE/////////////////////////
                                       
                                        for (int l = 0; l < FSMStates - 1; l++)
                                        {
                                            for (int m = l + 1; m < FSMStates; m++)
                                            {
                                                int input = getSepInfo(l, m);
                                                if (input >= 0)
                                                {
                                                    pair<int, int> st;
                                                    st.first = l;
                                                    st.second = m;
                                                    pair< pair<int, int>, vector<int>> sep;
                                                    sep.first = st;
                                                    vector<int> in;
                                                    in.push_back(input);
                                                    CS[0].push_back(in);
                                                    sep.second = in;
                                                    splittingInfo.push_back(sep);
                                                    int index = returnPairIndex(l, m, observedOutputs.size());
                                                    if (separatingSequences[index].size() == 0)
                                                    {
                                                        pairCount--;
                                                        separatingSequences[index].push_back(in);
                                                    }
                                                    for (int n = l; n < FSMStates - 1; n++)
                                                    {
                                                        for (int o = m + 1; o < FSMStates; o++)
                                                        {
                                                            if (separates(input, n, o))
                                                            {
                                                                pair<int, int> sts;
                                                                sts.first = n;
                                                                sts.second = o;
                                                                pair< pair<int, int>, vector<int>> seps;
                                                                seps.first = st;
                                                                seps.second = in;
                                                                splittingInfo.push_back(seps);
                                                                int indexs = returnPairIndex(n, o, observedOutputs.size());
                                                                if (separatingSequences[indexs].size() == 0)
                                                                {
                                                                    pairCount--;
                                                                    separatingSequences[indexs].push_back(in);
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        while (pairCount > 0)
                                        {
                                            for (int l = 0; l < FSMStates - 1; l++)
                                            {
                                                for (int m = l + 1; m < FSMStates; m++)
                                                {
                                                    int index = returnPairIndex(l, m, observedOutputs.size());
                                                    if (separatingSequences[index].size() == 0)
                                                    {
                                                        vector<int> divider;
                                                        for (int n = 0; n < FSMInputs; n++)
                                                        {
                                                            int ns1 = nextStates[l][n];
                                                            int ns2 = nextStates[m][n];
                                                            if (ns1 > ns2)
                                                            {
                                                                int t = ns1;
                                                                ns1 = ns2;
                                                                ns2 = t;
                                                            }

                                                            if (ns1 == ns2)
                                                            {
                                                                ;
                                                            }
                                                            else
                                                            {
                                                                int nIndex = returnPairIndex(ns1, ns2, nextStates.size());
                                                                if (separatingSequences[nIndex].size() > 0)
                                                                {
                                                                    divider.push_back(n);
                                                                }
                                                            }

                                                        }
                                                        if (divider.size() > 0)
                                                        {
                                                            int input = divider[rand() % divider.size()];
                                                            int ns1 = nextStates[l][input];
                                                            int ns2 = nextStates[m][input];
                                                            if (ns1 > ns2)
                                                            {
                                                                int t = ns1;
                                                                ns1 = ns2;
                                                                ns2 = t;
                                                            }
                                                            if (ns1 == ns2)
                                                            {
                                                                ;
                                                            }
                                                            else
                                                            {
                                                                int nIndex = returnPairIndex(ns1, ns2, nextStates.size());
                                                                vector<int>prevSeq;
                                                                prevSeq.push_back(input);
                                                                for (int o = 0; o < separatingSequences[nIndex][0].size(); o++)
                                                                {
                                                                    prevSeq.push_back(separatingSequences[nIndex][0][o]);
                                                                }
                                                                if (prevSeq.size() == 0)
                                                                    int ne = 0;
                                                                pairCount--;
                                                                separatingSequences[index].push_back(prevSeq);
                                                                CS[0].push_back(prevSeq);
                                                                pair<int, int> st;
                                                                st.first = l;
                                                                st.second = m;
                                                                pair< pair<int, int>, vector<int>> sep;
                                                                sep.first = st;
                                                                sep.second = prevSeq;
                                                                splittingInfo.push_back(sep);
                                                            }
                                                        }
                                                        else
                                                            ;

                                                        divider.clear();
                                                        divider.shrink_to_fit();
                                                    }
                                                }
                                            }
                                        }

                                        separatingSequences.clear();
                                        splittingInfo.clear();
                                        splittingInfo.shrink_to_fit();
                                        separatingSequences.shrink_to_fit();
                                        removePrefixesClassic(CS);
                                        clearCSClassic(CS);
                                        reachingSequences = calculateReachingSequneces();

                                        CSM = CS;
                                        tis = calculateTransitionIdentifiationSequences(reachingSequences, CS);
                                        sis = calculateStateIdentificationSequences(reachingSequences, CS);
                                        int MtisS = 0;
                                        for (int fg = 0; fg < tis.size(); fg++)
                                            MtisS += tis[fg].size();
                                        int MsisS = 0;
                                        for (int fg = 0; fg < sis.size(); fg++)
                                            MsisS += sis[fg].size();
                                        int Mtiscount = tis.size();
                                        int Msiscount = sis.size();
                                        classes.clear();
                                        classes.shrink_to_fit();
                                        Origin.clear();
                                        Origin.shrink_to_fit();
                                        WM = calculateCharacterisingSet(sis, tis);
                                        double memoryM = getUsedMemoryMB();
                                        auto end_min_T_time = std::chrono::high_resolution_clock::now();
                                        auto classicTime_T3 = end_min_T_time - start_min_T_time;
                                        float t_T3 = classicTime_T3 / std::chrono::milliseconds(1);
                                        /////////////////////////////Min Ends Here
                                        int WOlength = getLength(WO);
                                        int WNlength = getLength(WN);
                                        int CSOlength = getLength(CSO[0]);
                                        int CSNlength = getLength(CSN[0]);
                                        int WMlength = getLength(WM);
                                        int CSMlength = getLength(CSM[0]);
                                        
                                      {

                                            writeStatFile << FSMid << " " << FSMStates << " " << FSMInputs << " " << FSMOutputs << " " << FSMSize << " " << t_T1 << " " << memoryO <<" "  << transferLengthO << " " <<1<< " " << OsisS << " " << otiscount + osiscount << " " << OsisS + OtisS << " OWA" << endl;
                                            writeStatFile << FSMid << " " << FSMStates << " " << FSMInputs << " " << FSMOutputs << " " << FSMSize << " " << t_T2 << " " << memoryC << " " << transferLengthC << " " << nsiscount << " " << NsisS << " " << ntiscount + nsiscount << " " << NtisS + NsisS << " CWA" << endl;
                                            writeStatFile << FSMid << " " << FSMStates << " " << FSMInputs << " " << FSMOutputs << " " << FSMSize << " " << t_T3 << " " << memoryM << " " << transferLengthC << " " << Msiscount << " " << MsisS << " " << Mtiscount + Msiscount << " " << MtisS + MsisS << " MWA" << endl;

                                            writeSeqFile << FSMid << " " << FSMStates << " " << FSMInputs << " " << FSMOutputs << " " << FSMSize << " OrderedW={";
                                            for (int m = 0; m < WO.size(); m++)
                                            {
                                                writeSeqFile << "{";
                                                for (int n = 0; n < WO[m].size(); n++)
                                                {
                                                    writeSeqFile << char(WO[m][n] + 97);
                                                }
                                                writeSeqFile << "}";
                                            }
                                            writeSeqFile << "} --- ";
                                            writeSeqFile << " ClassicW={";
                                            for (int m = 0; m < WN.size(); m++)
                                            {
                                                writeSeqFile << "{";
                                                for (int n = 0; n < WN[m].size(); n++)
                                                {
                                                    writeSeqFile << char(WN[m][n] + 97);
                                                }
                                                writeSeqFile << "}";
                                            }
                                            writeSeqFile << "} --- ";
                                            writeSeqFile << " MinW={";
                                            for (int m = 0; m < WM.size(); m++)
                                            {
                                                writeSeqFile << "{";
                                                for (int n = 0; n < WM[m].size(); n++)
                                                {
                                                    writeSeqFile << char(WM[m][n] + 97);
                                                }
                                                writeSeqFile << "}";
                                            }
                                            writeSeqFile << " OCS={";
                                            for (int m = 0; m < CSO[0].size(); m++)
                                            {
                                                writeSeqFile << "{";
                                                for (int n = 0; n < CSO[0][m].size(); n++)
                                                {
                                                    writeSeqFile << char(CSO[0][m][n] + 97);
                                                }
                                                writeSeqFile << "}";
                                            }
                                            writeSeqFile << "} --- ";
                                            writeSeqFile << " CCS={";
                                            for (int m = 0; m < CSN[0].size(); m++)
                                            {
                                                writeSeqFile << "{";
                                                for (int n = 0; n < CSN[0][m].size(); n++)
                                                {
                                                    writeSeqFile << char(CSN[0][m][n] + 97);
                                                }
                                                writeSeqFile << "}";
                                            }
                                           writeSeqFile << "} --- ";
                                            writeSeqFile << " MCS={";
                                            for (int m = 0; m < CSM[0].size(); m++)
                                            {
                                                writeSeqFile << "{";
                                                for (int n = 0; n < CSM[0][m].size(); n++)
                                                {
                                                    writeSeqFile << char(CSM[0][m][n] + 97);
                                                }
                                                writeSeqFile << "}";
                                            }
                                            writeSeqFile << "}" << endl;
                                            counter++;
                                        }
                                        
                                        
                                        WO.clear();
                                        WO.shrink_to_fit();
                                        WN.clear();
                                        WN.shrink_to_fit();
                                    }
                                    CSO.clear();
                                    CSN.clear();
                                    CSO.shrink_to_fit();
                                    CSN.shrink_to_fit();
                                }
                                
                                CS.clear();
                                CS.shrink_to_fit();
                                activeStates.clear();
                                activeStates.shrink_to_fit();
                                indexes.clear();
                                indexes.shrink_to_fit();
                                edges.clear();
                                edges.shrink_to_fit();
                                oSIS.clear();
                                oSIS.shrink_to_fit();
                                
                                mySequence = stack<int>();
                                sequence = stack<int>();
                                myQ = queue<pair<int, vector<int>>>();
                                sms = vector < vector<bool>>();
                                indexes = deque<int>();
                                inputsOrder= vector<int>();
                                
                                edges= vector<vector<edge>>();
                                myQ = queue<pair<int, vector<int>>>();
                                sms = vector<vector<bool>>();;
                            }
                        
                    }   
                    writeSeqFile.close();
                    writeStatFile.close();    
}
