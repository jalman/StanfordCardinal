// FFT
const int g=137,bigp=15*(1<<27)+1,x=27;
vector<int> a,b;char ch;int C,i,N,L,len;
int powc(int a,int b){int d;if (b==0) return 1;d=powc(a,b/2);d=(LL)d*d%bigp;if (b%2==1) d=(LL)d*a%bigp;return d;}
vector<int> reverse(vector<int> a){
    int s,e,N=a.size();
    for (int i=0;i<N;i++){
        s=i,e=0;
        for (int j=1;j<=L;j++) e=(e<<1)+(s&1),s=s>>1;
        if (e<i) swap(a[e],a[i]);
    }
    return a;
}
vector<int> calc(vector<int> a,int o){
    int w,q,wi,x,y;
    for (int k=1;k<=L;k++){
        wi=powc(o,N>>k);
        for (int i=0;i<(N>>k);i++){
            w=1;
            for (int j=0;j<(1<<(k-1));j++){
                q=i*(1<<k)+j,x=a[q],y=(LL)(a[q+(1<<(k-1))])*w%bigp;
                a[q]=((LL)x+y)%bigp,a[q+(1<<(k-1))]=((LL)x-y+bigp)%bigp;
                w=(LL)w*wi%bigp;
            }
        }
    }
    return a;
}
vector<int> dft(vector<int> a){
    a=reverse(a);return calc(a,powc(g,(1<<x)/N));
}
vector<int> idft(vector<int> a){
    int m;a=reverse(a);a=calc(a,powc(g,(1<<x)/N*(N-1)));m=powc(N,bigp-2);
    for (int i=0;i<N;i++) a[i]=(LL)a[i]*m%bigp;
    return a;
}
int main()
{
    scanf("%d",&C);scanf("%c",&ch);a.clear();b.clear();
    for (int i=1;i<=C;i++){a.p_b(0);b.p_b(0);}
    for (int i=1;i<=C;i++){scanf("%c",&ch);a[C-i]=ch-48;}scanf("%c",&ch);
    for (int i=1;i<=C;i++){scanf("%c",&ch);b[C-i]=ch-48;}
    L=0,N=1;while ((N>>1)<C) N*=2,++L;
    while (a.size()<N) a.p_b(0);while (b.size()<N) b.p_b(0);
    /*core*/a=dft(a);b=dft(b);for (int i=0;i<N;i++) a[i]=(LL)a[i]*b[i]%bigp;a=idft(a);
    for (int i=0;i<N;i++){a[i+1]=(a[i+1]+a[i]/10);a[i]%=10;}
    len=N-1;while ((len>1)&&(a[len]==0)) --len;
    for (int i=len;i>=0;i--) printf("%d",a[i]);printf("\n");
}
