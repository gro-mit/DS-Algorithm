class GrayCode{
public:
    void graycode(int n,vector<string>&gray){
        if(n==1){
            gray.push_back("0");
            gray.push_back("1");
            return;
        }
        vector<string> lastgray=getGray(n-1);
        for(int i=0;i<lastgray.size();i++){
            gray.push_back("0"+lastgray[i]);
        }
        for(int i=lastgray.size()-1;i>=0;i--){
            gray.push_back("1"+lastgray[i]);
        }
        return;
    }
    vector<string> getGray(int n){
        vector<string> gray;
        graycode(n,gray);
        return gray;
    }
};