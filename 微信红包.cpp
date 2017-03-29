class Gift {
public:
    int getValue(vector<int> gifts, int n) {
        // write code here
        if(gifts.size()!=n||n==0)return 0;
        unordered_map<int,int> mp;
       //利用哈希表遍历计数，随着遍历的过程总会出现某红包值大于一半的时刻，则返回该红包。
        for(int i=0;i<n;i++){
            mp[gifts[i]]++;
            if(mp[gifts[i]]>n/2) return gifts[i];
        }
        //否则返回0
        return 0;
    }
};