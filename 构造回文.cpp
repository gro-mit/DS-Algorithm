/*
给定一个字符串s，你可以从中删除一些字符，使得剩下的串是一个回文串。如何删除才能使得回文串最长呢？
输出需要删除的字符个数。
*/

#include <bits/stdc++.h>
using namespace std;
#define MAX 101
int LCS(string& s,string& t){
    int ls=s.length(),lt=t.length();
    int C[MAX][MAX];
    for(int i=0;i<ls;i++){
        C[i][0]=0;
    }
    for(int i=0;i<lt;i++){
        C[0][i]=0;
    }
    for(int i=1;i<=ls;i++){
        for(int j=1;j<=lt;j++){
            if(s[i-1]==t[j-1]){
                C[i][j]=C[i-1][j-1]+1;
            } else{
                C[i][j]=max(C[i-1][j],C[i][j-1]);
            }
        }
    }
    return C[ls][lt];
}

int main(){
    string s;
    while(cin>>s){
        if(s.length()==1){
            cout<<1<<endl;
            continue;
        }
        string t=s;
        reverse(t.begin(),t.end());
        int res=s.length()-LCS(s,t);
        cout<<res<<endl;

    }
    return 0;
}