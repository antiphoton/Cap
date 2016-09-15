const double Z_INTERVAL = 0.5;
const double Z_HALFDEPTH = Z_INTERVAL/2;
const double ANNULUS_AREA = 20;
const double CIRCLE_BOTTOM = 8;

#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include<algorithm>
#include<fstream>
#include<map>
#include<string>
#include<vector>
using namespace std;

namespace ArgParser {
    map<string,string> m;
    void init(int argc,char **argv) {
        for (int i=0;i<argc;i++) {
            char *s=argv[i];
            if (strlen(s)>=3&&strncmp(s,"--",2)==0&&s[2]!='-') {
                if (i+1<argc&&argv[i+1][0]!='-') {
                    m[s+2]=argv[i+1];
                }
                else {
                    m[s+2]="";
                }
            }
        }
    };
    const char *const getString(const char *const name,const char *const defaultValue=0) {
        map<string,string>::iterator v=m.find(name);
        if (v==m.end()) {
            assert(defaultValue!=0);
            return defaultValue;
        }
        return v->second.c_str();
    }
    bool getBool(const string &name) {
        return m.find(name)!=m.end();
    }
    double getDouble0(const char *const name,double *p) {
        map<string,string>::iterator v=m.find(name);
        if (v==m.end()) {
            assert(p!=0);
            return *p;
        }
        const char *s=v->second.c_str();
        int n;
        double y;
        sscanf(s,"%lf%n",&y,&n);
        assert(n==strlen(s));
        return y;
    }
    double getDouble(const char *const name) {
        return getDouble0(name,0);
    }
    double getDouble(const char *const name,double defaultValue) {
        return getDouble0(name,&defaultValue);
    }
}
namespace Boundary {
    double xLow,xHigh,yLow,yHigh,zLow,zHigh;
    bool parseData0(const char **const s,const char *const sL,const char *const sH,double *pL,double *pH) {
        if (strcmp(s[2],sL)==0&&strcmp(s[3],sH)==0) {
            int n;
            sscanf(s[0],"%lf%n",pL,&n);
            bool bL=n==strlen(s[0]);
            sscanf(s[1],"%lf%n",pH,&n);
            bool bH=n==strlen(s[1]);
            return bL&&bH;
        }
        return false;
    }
    bool parseData(const char *const filename) {
        if (strlen(filename)==0) {
            return false;
        }
        FILE *f=fopen(filename,"r");
        static char line[100000];
        bool gX=false,gY=false,gZ=false;
        while (fgets(line,sizeof(line),f)) {
            static char s[4][1000];
            const char *p[]={s[0],s[1],s[2],s[3]};
            if (sscanf(line,"%s%s%s%s",s[0],s[1],s[2],s[3])!=4) {
                continue;
            }
            gX=gX||parseData0(p,"xlo","xhi",&xLow,&xHigh);
            gY=gY||parseData0(p,"ylo","yhi",&yLow,&yHigh);
            gZ=gZ||parseData0(p,"zlo","zhi",&zLow,&zHigh);
        }
        fclose(f);
        printf("%f\t%f\n",xLow,xHigh);
        printf("%f\t%f\n",yLow,yHigh);
        printf("%f\t%f\n",zLow,zHigh);
        return gX&&gY&&gZ;
    }
}
const double PI=acos(-1.0);
const double DENSITY_IDEAL = 0.033428;
double getRadius(double i) {
    return sqrt(ANNULUS_AREA*i/PI);
}
struct Position {
    double x,y,z;
    void read(FILE *f) {
        fscanf(f,"%lf%lf%lf",&x,&y,&z);
    };
    double distXY() const {
        return sqrt(x*x+y*y);
    };
};
class Snapshot {
public:
    int n;
    Position *p;
    Snapshot(FILE *f,const char *const waterType):p(0) {
        int n0;
        if (fscanf(f,"%d%*c",&n0)==EOF) {
            return;
        }
        Position *p0=new Position[n0];
        bool *valid=new bool[n0];
        fscanf(f,"%*[^\n]%*c");
        n=0;
        int i;
        for (i=0;i<n0;i++) {
            static char curType[10];
            fscanf(f,"%s",curType);
            valid[i]=strcmp(curType,waterType)==0;
            p0[i].read(f);
            if (valid[i]) {
                n++;
            }
        }
        p=new Position[n];
        n=0;
        for (i=0;i<n0;i++) {
            if (valid[i]) {
                memcpy(&p[n],&p0[i],sizeof(Position));
                n++;
            }
        }
        delete[] valid;
        delete[] p0;
        if (1) {
            handleBoundary((int)(long)(char*)(&(((Position*)0)->x)),Boundary::xLow,Boundary::xHigh);
            handleBoundary((int)(long)(char*)(&(((Position*)0)->y)),Boundary::yLow,Boundary::yHigh);
        }
        handleMassCenter((int)(long)(char*)(&(((Position*)0)->x)));
        handleMassCenter((int)(long)(char*)(&(((Position*)0)->y)));
    }
    ~Snapshot() {
        delete[] p;
    };
private:
    double getNearest(double w,int b) {
        double ret;
        int i;
        for (i=0;i<n;i++) {
            double cur=fabs(w-*(double*)((char*)&p[i]+b));
            if (i==0||cur<ret) {
                ret=cur;
            }
        }
        return ret;
    };
    void handleBoundary(int b,double low,double high) {
        double max=-1;
        double best;
        double w;
        double width=high-low;
        for (w=low;w<high;w+=width/20) {
            double cur=getNearest(w,b);
            if (cur>max) {
                max=cur;
                best=w;
            }
        }
        double center=best+width/2;
        if (center>=high) {
            center-=width;
        }
        for (int i=0;i<n;i++) {
            double *px=(double *)((char *)&p[i]+b);
            *px-=center;
            if (*px>=width/2) {
                *px-=width;
            }
        }
    };
    void handleMassCenter(int b) {
        double sum=0;
        for (int i=0;i<n;i++) {
            double *px=(double *)((char *)&p[i]+b);
            sum+=*px;
        }
        double average=sum/n;
        for (int i=0;i<n;i++) {
            double *px=(double *)((char *)&p[i]+b);
            *px-=average;
        }
    };
};
inline double sqr(double x) {
    return x*x;
}
pair<double,double> fitCircle(const vector< pair<double,double> > &v) {
    int n=v.size();
    double h,r;
    {
        double x1=v[0].first,y1=v[0].second;
        double y2=v[n-1].second;
        h=(sqr(x1)+sqr(y1)-sqr(y2))/(2*(y1-y2));
        r=(sqr(x1)+sqr(y2-y1))/(2*(y2-y1));
    }
    int it=10;
    while (it--) {
        double h11=0,h12=0,h22=0;
        double g1=0,g2=0;
        for (int i=0;i<n;i++) {
            double x=v[i].first;
            double y=v[i].second;
            double d=sqrt(sqr(x)+sqr(y-h));
            h11+=2-2*r*sqr(x)/(d*d*d);
            h12+=2*(y-h)/d;
            h22+=2;
            g1+=-2*(y-h)*(d-r)/d;
            g2+=-2*(d-r);
        }
        double deltaH=h11*h22-sqr(h12);
        double ht11=h22/deltaH;
        double ht12=-h12/deltaH;
        double ht22=h11/deltaH;
        double m1=g1*ht11+g2*ht12;
        double m2=g1*ht12+g2*ht22;
        h-=m1;
        r-=m2;
    }
    return make_pair(h,r);
}
vector<Snapshot *> trajectory;
int main(int argc,char **argv) {
    assert(1/(1+sqr(1.0/0))==0);
    ArgParser::init(argc,argv);
    assert(
        Boundary::parseData(ArgParser::getString("boundary-data",""))
    );
    const char *const title=ArgParser::getString("title","");
    const double pBegin=ArgParser::getDouble("frame-begin",0.5);
    FILE *fGnuplot=fopen("./plot.gp","w");
    FILE *fXyz=fopen(ArgParser::getString("xyz"),"r");
    const char *const waterType=ArgParser::getString("water");
    const int ATOMS_PER_NOTIFICATION=1000000;
    long long atomParsed=0;
    while (1) {
        Snapshot *s=new Snapshot(fXyz,waterType);
        if (s->p==0) {
            delete s;
            break;
        }
        trajectory.push_back(s);
        atomParsed+=s->n;
        if (atomParsed>=ATOMS_PER_NOTIFICATION) {
            printf("Snapshots parsed: %d\n",trajectory.size());
            atomParsed%=ATOMS_PER_NOTIFICATION;
        }
    }
    const int nFrame=trajectory.size();
    FILE *fSizeHistory=fopen("shapeHistory.txt","w");
    FILE *fEffectiveRadius=fopen("effectiveRadius.txt","w");
    const int iBegin=(int)(pBegin*nFrame);
    const int iEnd=nFrame-1;
    for (int i=0;i<nFrame;i++) {
        double z=0,r=0,sC=0,sS=0;
        const Position *p=trajectory[i]->p;
        int n=trajectory[i]->n;
        for (int j=0;j<n;j++) {
            z+=p[j].z;
            r+=sqr(p[j].x)+sqr(p[j].y);
            double t=atan2(p[j].y,p[j].x);
            sC+=cos(4*t);
            sS+=sin(4*t);
        }
        z/=n;
        r/=n;
        r=sqrt(r);
        sC/=n;
        sS/=n;
        fprintf(fSizeHistory,"%f\t%f\t%f\n",z,r,sqrt(sqr(sC)+sqr(sS)));
    }
    fprintf(fGnuplot,"set terminal svg size 1600,900;\n");
    fprintf(fGnuplot,"set output \'shapeHistory.svg\';\n");
    fprintf(fGnuplot,"set xrange [%d:%d];\n",-1,nFrame);
    fprintf(fGnuplot,"set xtics rotate;\n",-1,nFrame);
    fprintf(fGnuplot,"set multiplot layout 1,3 title \'%s\';\n",title);
    fprintf(fGnuplot,"set title \'Average height\';\n");
    fprintf(fGnuplot,"plot \'shapeHistory.txt\' ");
    fprintf(fGnuplot,"u 0:(($0<%d||$0>%d)?$1:(1/0)) t \'\' w l lt rgb \'#ff0000\', ",iBegin,iEnd);
    fprintf(fGnuplot,"\'\' u 0:(($0>=%d&&$0<%d)?$1:(1/0)) t \'\' w l lt rgb \'#0000ff\';\n",iBegin,iEnd);
    fprintf(fGnuplot,"set title \'Average radius\';\n");
    fprintf(fGnuplot,"plot \'shapeHistory.txt\' ");
    fprintf(fGnuplot,"u 0:(($0<%d||$0>%d)?$2:(1/0)) t \'\' w l lt rgb \'#ff0000\', ",iBegin,iEnd);
    fprintf(fGnuplot,"\'\' u 0:(($0>=%d&&$0<%d)?$2:(1/0)) t \'\' w l lt rgb \'#0000ff\';\n",iBegin,iEnd);
    fprintf(fGnuplot,"set title \'Order\';\n");
    fprintf(fGnuplot,"plot \'shapeHistory.txt\' ");
    fprintf(fGnuplot,"u 0:(($0<%d||$0>%d)?($3/%f):(1/0)) t \'\' w l lt rgb \'#ff0000\', ",iBegin,iEnd,PI-3);
    fprintf(fGnuplot,"\'\' u 0:(($0>=%d&&$0<%d)?($3/%f):(1/0)) t \'\' w l lt rgb \'#0000ff\';\n",iBegin,iEnd,PI-3);
    fprintf(fGnuplot,"unset multiplot;\n");
    fprintf(fGnuplot,"set xtics norotate;\n",-1,nFrame);
    vector< pair<double,double> > effectiveRadiusList;
    for (double zCenter=Boundary::zLow;zCenter<Boundary::zHigh;zCenter+=Z_INTERVAL) {
        printf("zCenter=%f\n",zCenter);
        double z1=zCenter-Z_HALFDEPTH,z2=zCenter+Z_HALFDEPTH;
        vector<double> v;
        for (int i=iBegin;i<=iEnd;i++) {
            const Position *p=trajectory[i]->p;
            for (int j=trajectory[i]->n-1;j>=0;j--) {
                if (p[j].z>z1&&p[j].z<z2) {
                    v.push_back(p[j].distXY());
                }
            }
        }
        sort(v.begin(),v.end());
        double effectiveRadius=-1;
        bool solid0=true,solid2;
        int j=0;
        for (int i=1;;i++) {
            const double rCur=getRadius(i);
            int c=0;
            while (1) {
                if (j>=v.size()||v[j]>rCur) {
                    break;
                }
                c++;
                j++;
            }
            double d=1.0*c/(iEnd-iBegin+1)/(ANNULUS_AREA*Z_HALFDEPTH*2)/DENSITY_IDEAL;
            solid2=d>0.5;
            if (solid0==true&&solid2==false) {
                effectiveRadius=getRadius(i-1);
            }
            if (c==0&&effectiveRadius>=0) {
                if (solid0==false) {
                    fprintf(fEffectiveRadius,"%f\t%f\n",zCenter,effectiveRadius);
                    if (zCenter>=CIRCLE_BOTTOM) {
                        if (effectiveRadiusList.size()==0||effectiveRadiusList.back().first>0) {
                            effectiveRadiusList.push_back(make_pair(effectiveRadius,zCenter));
                        }
                    }
                    break;
                }
            }
            solid0=solid2;
        }
    }
    fprintf(fGnuplot,"set terminal svg size 1600,900;\n");
    fprintf(fGnuplot,"set output \'effectiveRadius.svg\';\n");
    fprintf(fGnuplot,"set title \'%s';\n",title);
    fprintf(fGnuplot,"set xrange [0:%f];\n",(Boundary::yHigh-Boundary::yLow)/2);
    fprintf(fGnuplot,"set yrange [%f:%f];\n",Boundary::zLow,Boundary::zHigh);
    fprintf(fGnuplot,"set size ratio -1;\n");
    fprintf(fGnuplot,"plot \'effectiveRadius.txt\' u 2:1 t \'\',");
    {
        pair<double,double> circleParam=fitCircle(effectiveRadiusList);
        double h=circleParam.first;
        double r=circleParam.second;
        fprintf(fGnuplot,"sqrt(%f*%f-x*x)%s%f t \'cos(theta)=%f\';\n",r,r,h>=0?"+":"",h,-h/r);
    }
    for (int i=0;i<trajectory.size();i++) {
        delete trajectory[i];
    }
    return 0;
}

