#include <iostream>
#include <iomanip>
#include <algorithm>
#include <time.h>
#include <random>

#define SIDE_LENGTH 100
#define EXPTIME 5
#define EPS 1e-6
using namespace std;

struct Point{
    double x,y;
    bool status = false;
    Point():x(0),y(0){}
    Point(double ax, double ay):x(ax),y(ay){}
    bool operator < (const Point &p)
    {
        if(x-p.x < EPS)
        {
            return y < p.y;
        }
        return x < p.x;
    }
    bool operator == (const Point &p)
    {
        return (abs(x-p.x) <= EPS && abs(y-p.y) <= EPS);
    }
    bool operator != (const Point &p)
    {
        return (abs(x-p.x) > EPS || abs(y-p.y) > EPS);
    }
    void printPoint()
    {
        cout<<setiosflags(ios::fixed)<<setprecision(2);
        cout<<"->"<<'('<<x<<','<<y<<')';
    }
};

struct ExpResult{
    int method = 1;
    string methodName;
    int pointsNumber = 0;
    double runtime = 0;
    bool valid = false;
    ExpResult(int hullMethod, int num, double runningTime):method(hullMethod),pointsNumber(num),runtime(runningTime)
    {
        switch(hullMethod)
        {
            case 1:
                methodName = "BruteForce";
                break;
            case 2:
                methodName = "Graham-Scan";
                break;
            case 3:
                methodName = "Divide&Conquer";
                break;
            default:
                methodName = "";
        }
    }
};

//generate a set of random points
vector<Point> generateSampleSet(int n)
{
    default_random_engine xGenerator((unsigned int)time(NULL));
    default_random_engine yGenerator((unsigned int)time(NULL)+1);
    uniform_real_distribution<double> uniformDistribution(0.0,SIDE_LENGTH);
    double x,y;
    vector<Point> sampleSet;
    for(int i=0;i<n;++i){
        x = uniformDistribution(xGenerator);
        y = uniformDistribution(yGenerator);
        Point point(x,y);
        sampleSet.push_back(point);
    }
    return sampleSet;
}

//get the origin point
int getOriginPoint(vector<Point> &set)
{
    int index = 0;
    for(int i = 1;i<set.size();++i)
    {
        if(set[i] < set[index])
            index = i;
    }
    return index;
}

//calculate cross product
double crossProduct(Point a, Point b, Point p)
{
    //calculate cross product
    return (b.x-a.x)*(p.y-a.y)-(p.x-a.x)*(b.y-a.y);
}

//calculate the distance between 2 points
double distance(Point a, Point b)
{
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}

//judge which side of the line does the point belong to,if exists the same angle,select the point by their distance
bool cmp(Point &origin, Point &a, Point &b)
{
    double product = crossProduct(origin, a, b);
    if(abs(product) < EPS)
    {
        return distance(origin, a) < distance(origin, b);
    } else{
        //if >0, b is on the left side of line a-origin
        return product>0;
    }
}

//sort
int partition(vector<Point> &set, int l, int r)
{
    Point key = set[l];
    while(l < r)
    {
        while(l < r && cmp(set[0], key, set[r])) --r;
        swap(set[l], set[r]);
        while(l < r && cmp(set[0], set[l], key)) ++l;
        swap(set[l], set[r]);
    }
    return l;
}
void sort(vector<Point> &set, int l, int r)
{
    if(l < r)
    {
        int index = partition(set, l, r);
        sort(set, l, index-1);
        sort(set, index+1, r);
    }
}
////get the boundary of hull
//void getHullBoundary(vector<Point> &set, vector<Point> &res)
//{
//    for(int i = 0;i<set.size();++i)
//    {
//        if(set[i].status)
//        {
//            res.push_back(set[i]);
//        }
//    }
//}

//show the boundary of hull
void showHullBoundary(vector<Point> &set)
{
    int cnt = 0;
    for(int i = 0;i < set.size();++i)
    {
        if(set[i].status)
        {
            set[i].printPoint();
            ++cnt;
            if(cnt%5 == 0)
                cout<<endl;
        }
    }
    cout<<endl;
}

//judge the results of those three methods
bool judgeHullBoundaryResult(vector<Point> &set1, vector<Point> &set2, vector<Point> &set3, int n)
{
    for(int i = 0;i < n, set1[i].status&&set2[i].status&&set3[i].status; ++i)
    {
        if(set1[i] != set2[i] || set1[i] != set3[i] || set2[i] != set3[i]){
            cout<<"ERROR!"<<endl;
            return false;
        }
    }
    cout<<"THE SAME BOUNDARY."<<endl;
    return true;
}
//BruteForce way to solve convex hull
ExpResult bruteForceConvexHull(vector<Point> &set, int n)
{
    clock_t start,end;
    start = clock();
    for(int i = 0;i < n-1;++i)
    {
        for(int j = i+1;j < n;++j)
        {
            int cnt0 = 0,cnt1 = 0;
            for(int k = 0;k < n;++k)
            {
                if(k != i && k != j)
                {
                    if(crossProduct(set[i], set[j], set[k]) >= 0)
                    {
                        ++cnt0;
                    }
                    if(crossProduct(set[i],set[j],set[k]) <= 0)
                    {
                        --cnt1;
                    }
                }
            }
            if(cnt0 == n-2||cnt1 == 2-n)
            {
                set[i].status = true;
                set[j].status = true;
            }
        }
    }
    end = clock();
    double duration = (double)(end - start)*1000/CLOCKS_PER_SEC;
    sort(set, 1, n-1);
    showHullBoundary(set);
    cout<<n<<" points' runtime: "<<duration<<" ms"<<endl;
    return ExpResult(1, n, duration);
}

//Graham-Scan way to solve convex hull
ExpResult grahamScanConvexHull(vector<Point> &set, int n)
{
    clock_t start,end;
    start = clock();
    int top = 2;
    vector<Point> boundaryPoints;
    boundaryPoints.push_back(set[0]);
    sort(set,1,n-1);
    boundaryPoints.push_back(set[1]);
    boundaryPoints.push_back(set[2]);
    for(int i = 3;i < n;++i)
    {
        while(cmp(set[i],boundaryPoints[top],boundaryPoints[top-1]))
        {
            --top;
            boundaryPoints.pop_back();
        }
        boundaryPoints.push_back(set[i]);
        ++top;
    }
    for(int i = 0;i < n;++i)
    {
        for(int j = 0;j < boundaryPoints.size();++j)
        {
            if(set[i] == boundaryPoints[j])
                set[i].status = true;
        }
    }
    end = clock();
    double duration = (double)(end - start)*1000/CLOCKS_PER_SEC;
    showHullBoundary(set);
    cout<<n<<" points' runtime: "<<duration<<" ms"<<endl;
    return ExpResult(2, n, duration);
}

//Divide & Conquer way to solve convex hull
void divideConquer(vector<Point> &set, int start, int end, bool isLeft)
{
    if(start >= end) return;
    double maxDis = 0,minDis = 0;
    int maxIndex = -1, minIndex = -1;
    double distance = 0;
    for(int i = start+1;i < end;++i)
    {
        distance = crossProduct(set[start], set[end], set[i]);
        if(distance < minDis && !isLeft)
        {
            minDis = distance;
            minIndex = i;
        }
        if(distance > maxDis && isLeft)
        {
            maxDis = distance;
            maxIndex = i;
        }
    }
    if(minIndex > 0 && !isLeft)
    {
        set[minIndex].status = true;
        divideConquer(set, start, minIndex, isLeft);
        divideConquer(set, minIndex, end, isLeft);
    }
    if(maxIndex > 0 && isLeft)
    {
        set[maxIndex].status = true;
        divideConquer(set, start, maxIndex, isLeft);
        divideConquer(set, maxIndex, end, isLeft);
    }
}
ExpResult divideConquerConvexHull(vector<Point> &set, int n)
{
    clock_t start,end;
    start = clock();
    sort(set, 1, n-1);
    set[0].status = true;
    set[n-1].status = true;
    divideConquer(set, 0, n-1, true);
    divideConquer(set, 0, n-1, false);
    end = clock();
    double duration = (double)(end - start)*1000/CLOCKS_PER_SEC;
    showHullBoundary(set);
    cout<<n<<" points' runtime: "<<duration<<" ms"<<endl;
    return ExpResult(3, n, duration);
}

//experiment flow
void expFlow(vector<ExpResult> &res, int (&pointNumber)[EXPTIME], int n)
{
    for(int round = 0;round < n;++round)
    {
        vector<Point> set=generateSampleSet(pointNumber[round]),set1,set2,set3;
        int origin = getOriginPoint(set);
        swap(set[0], set[origin]);
        for(int i = 0;i < pointNumber[round];++i)
        {
            set1.push_back(set[i]);
            set2.push_back(set[i]);
            set3.push_back(set[i]);
        }
        cout<<endl;
        cout << "convex hull: bruteforce way" << endl;
        res.push_back(bruteForceConvexHull(set1, pointNumber[round]));
        cout<<"*********************************************"<<endl;
        cout << "convex hull: Graham-Scan way" << endl;
        res.push_back(grahamScanConvexHull(set2, pointNumber[round]));
        cout<<"*********************************************"<<endl;
        cout << "convex hull: Divide & Conquer way" << endl;
        res.push_back(divideConquerConvexHull(set3, pointNumber[round]));
        cout<<"*********************************************"<<endl;
        if(judgeHullBoundaryResult(set1, set2, set3, pointNumber[round]))
        {
            for(int i = 0;i < 3;++i)
            {
                res[3*round+i].valid = true;
            }
        }
    }
    for(int i = 0; i < n;++i)
    {
        cout<<pointNumber[i]<<" points:"<<endl;
        for(int j = 0;j < 3;++j)
        {
            if(res[i*3+j].valid){
                cout<<setiosflags(ios::fixed)<<setprecision(2);
                cout<<res[i*3+j].methodName<<": runtime is "<<res[i].runtime<<" ms"<<endl;
            } else{
                string excep = res[i*3+j].methodName+" ERROR!";
                throw excep;
            }
        }
        cout<<"+++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    }
}

int main()
{
    int n = EXPTIME;
    int pointNumber[EXPTIME] = {10, 20, 30, 50, 80};
    vector<ExpResult> res;
    try{
        expFlow(res, pointNumber, n);
    }catch(const char* msg) {
        cerr<<msg<<endl;
    }
    return 0;
}