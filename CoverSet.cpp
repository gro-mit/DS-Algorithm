#include <iostream>
#include <iomanip>
#include <algorithm>
#include <time.h>
#include <vector>
#include <set>
#include <map>
#include <random>
#include <Eigen/Dense>
#include "SimplexSolver.h"
//#include "exception.h"
#define MAXSIZE 5000
#define SUBSETSIZE 20
#define EXPTIME 3
using namespace std;
using namespace Eigen;

struct Point{
    int pointId;
    double pointWeight = -1;
    bool status = false;
    Point():pointId(-1){}
    Point(int id):pointId(id){}
    Point(int id, double weight):pointId(id),pointWeight(weight){}
    bool operator == (const Point &p)
    {
        return pointId == p.pointId;
    }
    bool operator != (const Point &p)
    {
        return pointId != p.pointId;
    }
    void printPoint()
    {
        cout<<"->"<<"Point Id: "<<pointId;
    }
};

//generate a set with n elements
vector<Point> generateSampleSet(int n)
{
    vector<Point> set;
    default_random_engine generator((unsigned int)time(NULL));
    uniform_real_distribution<double> uniformDistribution(0.01,1);
    for(int i = 0;i < n;++i)
    {
        double weight = uniformDistribution(generator);
        Point element(i, weight);
        set.push_back(element);
    }
    return set;
}

//generate subsets
int countIndex(bool index[], int n)
{
    int cnt = 0;
    for(int i = 0;i < n;++i)
    {
        if(!index[i])
            ++cnt;
    }
//    cout<<cnt<<" not selected"<<endl;
    return cnt;
}

int countIndex(vector<Point> &set)
{
    int cnt = 0;
    for(int i = 0;i < set.size();++i)
    {
        if(!set[i].status)
            ++cnt;
    }
//    cout<<cnt<<" left"<<endl;
    return cnt;
}

//set to vector
void convertSetToVector(set<int> &sampleSet, vector<int> &sampleVector)
{
    vector<int> copyVector(sampleSet.begin(), sampleSet.end());
    sampleVector.clear();
    sampleVector = copyVector;
}

vector<vector<int>> generateSubSet(int n)
{
    default_random_engine generator((unsigned int)time(NULL));
    uniform_int_distribution<int> uniformDistribution(0,n-1);
    bool index[MAXSIZE] = {false};
    int id;
    vector<vector<int>> subsets;
    vector<int> subset;
    set<int> firstSet;
    while(firstSet.size() != SUBSETSIZE)
    {
        id = uniformDistribution(generator);
        index[id] = true;
        firstSet.insert(id);
    }
    convertSetToVector(firstSet, subset);
    firstSet.clear();
    subsets.push_back(subset);
    subset.clear();
    uniform_int_distribution<int> uniformDistributionSubset(0,SUBSETSIZE);//generate a random number in [0,20]
    int sizeOfNextSubset, newElementNumber, oldElementNumber, idInSubset;
    while(countIndex(index, n) >= SUBSETSIZE)
    {
        sizeOfNextSubset = uniformDistributionSubset(generator);
        if(sizeOfNextSubset == 0)
        {
//            cout<<"sizeOfNextSubset is 0, continue"<<endl;
            continue;
        }
        uniform_int_distribution<int> uniformDistributionNextSubset(0,sizeOfNextSubset);
        newElementNumber = uniformDistributionNextSubset(generator);
        oldElementNumber = sizeOfNextSubset - newElementNumber;
        if(newElementNumber == 0)
        {
//            cout<<"newElementNumber is 0, continue"<<endl;
            continue;
        }
        set<int> tempSet;
        while(tempSet.size() != newElementNumber)
        {
            idInSubset = uniformDistribution(generator);
            if(!index[idInSubset])
            {
                index[idInSubset] = true;
                tempSet.insert(idInSubset);
            }
        }
        while(tempSet.size() != sizeOfNextSubset)
        {
            vector<int> lastSet = subsets[subsets.size()-1];
            if(lastSet.size() <= oldElementNumber)
            {
                convertSetToVector(tempSet, subset);
                subsets.push_back(subset);
                subset.clear();
                break;
            }else{
                uniform_int_distribution<int> uniformDistributionLastSet(0, lastSet.size()-1);
                idInSubset = uniformDistributionLastSet(generator);
                idInSubset = lastSet[idInSubset];
                tempSet.insert(idInSubset);
            }
        }
        convertSetToVector(tempSet, subset);
        subsets.push_back(subset);
        subset.clear();
    }
    for(int i = 0;i < n;++i)
    {
        if(!index[i])
        {
            index[i] = true;
            subset.push_back(i);
        }
    }
    subsets.push_back(subset);
    subset.clear();
    int randomSubsetNumber = n/SUBSETSIZE;
    for(int i = 0;i < randomSubsetNumber;++i)
    {
        sizeOfNextSubset = uniformDistributionSubset(generator);//generate random subset size
        set<int> tempSet;
        while(tempSet.size() != sizeOfNextSubset)
        {
            id = uniformDistribution(generator);
            tempSet.insert(id);
        }
        convertSetToVector(tempSet, subset);
        subsets.push_back(subset);
        subset.clear();
    }
    return subsets;
}

//sample invalid testing
bool invalidTesting(vector<vector<int>> &samples, int n)
{
    vector<int>::iterator maxId,minId;
    for(int i = 0;i < samples.size();++i)
    {
        maxId = max_element(begin(samples[i]), end(samples[i]));
        minId = min_element(begin(samples[i]), end(samples[i]));
        if(*maxId > n || *minId < 0)
        {
            cout<<"SAMPLE INVALID TESTING: ERROR SAMPLE!"<<endl;
            return false;
        }
        cout<<"SAMPLE INVALID TESTING: ALL SAMPLE PASS."<<endl;
        return true;
    }
}

//calculate cover scale
int findMaxCoverScale(vector<Point> &sampleSet, vector<vector<int>> &subsets)
{
    int maxUncovered = -1, index = 0;
    for(int i = 0;i < subsets.size();++i)
    {
        int cnt = 0;
        for(int j = 0;j < subsets[i].size();++j)
        {
            if(!sampleSet[subsets[i][j]].status)
            {
                ++cnt;
            }
        }
        if(maxUncovered < cnt)
        {
            maxUncovered = cnt;
            index = i;
        }
    }
    return index;
}

//show cover
void showCover(vector<vector<int>> &sets, int n)
{
    set<int> res;
    for(int i = 0;i < sets.size();++i)
    {
        for(int j = 0;j < sets[i].size();++j)
        {
            res.insert(sets[i][j]);
        }
    }
    if(res.size() == n)
    {
        cout<<setiosflags(ios::fixed)<<setprecision(2);
        cout<<"ALL SETS COVERED"<<endl;
    } else{
        cout<<"ERROR RESULT!"<<endl;
    }
}

//greedy approx
void greedySetCover(vector<Point> &sampleSet, vector<vector<int>> &subsets, vector<vector<int>> &answer)
{
    clock_t start,end;
    start = clock();
    while(countIndex(sampleSet) != 0)
    {
        int maxCoverIndex = findMaxCoverScale(sampleSet, subsets);
        for(int i = 0;i < subsets[maxCoverIndex].size();++i)
        {
            sampleSet[subsets[maxCoverIndex][i]].status = true;
        }
        answer.push_back(subsets[maxCoverIndex]);
    }
    end = clock();
    double duration = (double)(end - start)*1000/CLOCKS_PER_SEC;
    showCover(answer, sampleSet.size());
    cout<<sampleSet.size()<<" points' runtime: "<<duration<<" ms"<<endl;
}

//find element with highest occurence frequency
int findHighestElement(vector<vector<int>> &set)
{
    map<int,int> mp;
    int highOcc = -1;
    for(int i = 0;i < set.size();++i)
    {
        for(int j = 0;j < set[i].size();++j)
        {
            mp[(set[i])[j]]++;
            if(mp[(set[i])[j]] > highOcc)
                highOcc = mp[(set[i])[j]];
        }
    }
    return highOcc;
}

//calculate weight of subset
VectorXd calculateSubsetWeight(vector<Point> &sampleSet, vector<vector<int>> &subsets)
{
    VectorXd weights(subsets.size());

    for(int i = 0;i < subsets.size();++i)
    {
        double weight = 0;
        for(int j = 0;j < subsets[i].size();++j)
        {
            weight += sampleSet[subsets[i][j]].pointWeight;
        }
        weight /= subsets[i].size();
        weights(i) = weight;
    }
    return weights;
}

//calculate coefficient matrix
MatrixXd calculateCoefficientMatrix(vector<Point> &sampleSet, vector<vector<int>> &subsets)
{
    MatrixXd A(subsets.size(), sampleSet.size());
    A.setZero();
    for(int i = 0;i < subsets.size();++i)
    {
        for(int j = 0;j < subsets[i].size();++j)
        {
            A(i,subsets[i][j]) = 1;
        }
    }
    return A;
}

//get coefficient c
VectorXd getCofficientC(int n)
{
    VectorXd res(n);
    res.setOnes();
    return res;
}

VectorXd getCofficientC(int n, int ones)
{
    VectorXd res(n);
    res.setOnes();
    for(int i = ones;i < n;++i)
    {
        res(i) = 0;
    }
    return res;
}

//linear programming
void linearProgrammingSetCover(vector<Point> &sampleSet, vector<vector<int>> &subsets, vector<vector<int>> &answer)
{
    clock_t start,end;
    start = clock();
    VectorXd weights = calculateSubsetWeight(sampleSet, subsets);
    VectorXd baseVector = getCofficientC(subsets.size());
    MatrixXd baseMatrix(baseVector.asDiagonal());
    MatrixXd coefficient = calculateCoefficientMatrix(sampleSet, subsets);
    MatrixXd coefficientMatrix(coefficient.rows(), coefficient.cols()+baseMatrix.cols());
    coefficientMatrix << coefficient, baseMatrix;
    VectorXd c = getCofficientC(coefficientMatrix.cols(), sampleSet.size());

    VectorXd objFunction(sampleSet.size());
    objFunction.setOnes();
    MatrixXd constraints(coefficient.rows(),coefficient.cols()+weights.cols());
    constraints << coefficient, weights;

    SimplexSolver* LPSolver = NULL;
    LPSolver = new SimplexSolver(SIMPLEX_MAXIMIZE, objFunction, constraints);
    if(LPSolver->hasSolution())
    {
        cout<<"MAX is: "<<LPSolver->getOptimum()<<endl;
        cout<<"Got ys' SOLUTION"<<endl;
    } else{
        cout<<"NO SOLUTION!"<<endl;
    }
    VectorXd solution = LPSolver->getSolution();
    while(countIndex(sampleSet) != 0)
    {
        for(int i = 0;i < sampleSet.size();++i){
            if(!sampleSet[i].status)
            {
                int whichSubset = 0;
                double minWeight = INFINITY;
                vector<int>::iterator findSample;
                for(int j = 0;j < subsets.size();++j)
                {
                    if(find(subsets[j].begin(), subsets[j].end(), sampleSet[i].pointId) != subsets[j].end())
                    {
                        if(weights(j) < minWeight)
                        {
                            minWeight = weights(j);
                            whichSubset = j;
                        }
                    }
                }
                for(int k = 0;k < subsets[whichSubset].size();++k)
                {
                    sampleSet[subsets[whichSubset][k]].status = true;
                }
                answer.push_back(subsets[whichSubset]);
            }
        }
    }
    end = clock();
    double duration = (double)(end - start)*1000/CLOCKS_PER_SEC;
    showCover(answer, sampleSet.size());
    cout<<sampleSet.size()<<" points' runtime: "<<duration<<" ms"<<endl;
}

//cmp two answers
void compareTwoAnswers(vector<vector<int>> &ans1, vector<vector<int>> &ans2)
{
    map<int,int> cmp;
    if(ans1.size() != ans2.size())
    {
        cout<<"SIZE NOT MATCHED, ERROR ANSWERS!"<<endl;
        return;
    }
    for(int i = 0;i < ans1.size();++i){
        ++cmp[ans1.size()];
        ++cmp[ans2.size()];
    }
    map<int,int>::iterator it = cmp.begin();
    while(it != cmp.end())
    {
        if((it++)->second != 2)
        {
            cout<<"ERROR ANSWERS!"<<endl;
            return;
        }
    }
    cout<<"ANSWERS RIGHT"<<endl;
}

//experiment flow
void expFlow(int (&pointNumber)[EXPTIME], int n)
{
    for(int round = 0;round < n;++round)
    {
        cout<<"+++++++++++++++++ ROUND "<<round+1<<" ++++++++++++++++++"<<endl;
        vector<Point> sampleSet = generateSampleSet(pointNumber[round]);
        vector<Point> set1 = sampleSet, set2 = sampleSet;
        vector<vector<int>> subsets = generateSubSet(pointNumber[round]);
        vector<vector<int>> subsets1 = subsets, subsets2 = subsets;
        vector<vector<int>> answer1, answer2;
        cout<<"set cover: greedy way"<<endl;
        greedySetCover(set1, subsets2, answer1);
        cout<<"********************************************"<<endl;
        cout<<"set cover: linear programming way"<<endl;
        linearProgrammingSetCover(set2, subsets2, answer2);
        cout<<"********************************************"<<endl;
    }
}

int main()
{
    int n = EXPTIME;
    int pointNumber[EXPTIME] = {100, 1000, 5000};
    try{
        expFlow(pointNumber, n);
    }catch(const char* msg){
        cerr<<msg<<endl;
    }
    return 0;
}