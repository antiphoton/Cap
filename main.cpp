const bool BOUNDARY = true;
const double X_LOW = 0;
const double X_HIGH = 199.22048388651729;
const double Y_LOW = 0;
const double Y_HIGH = 200.22;
const double Z_LOW = 0;
const double Z_HIGH = 62;
const char *WATER_TYPE = "1";
const double Z_INTERVAL = 0.5;
const double Z_HALFDEPTH = Z_INTERVAL/2;
const double ANNULUS_AREA = 20;
const double FRAME_BEGIN = 0.5;
const double CIRCLE_BOTTOM = 8;
const char *TITLE = "title";

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
    Snapshot(FILE *f):p(0) {
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
            valid[i]=strcmp(curType,WATER_TYPE)==0;
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
        if (BOUNDARY) {
            handleBoundary((int)(long)(char*)(&(((Position*)0)->x)),X_LOW,X_HIGH);
            handleBoundary((int)(long)(char*)(&(((Position*)0)->y)),Y_LOW,Y_HIGH);
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
    ArgParser::init(argc,argv);
    FILE *gnuplot=fopen("./plot.gp","w");
    FILE *fXyz=fopen(ArgParser::getString("xyz"),"r");
    const int ATOMS_PER_NOTIFICATION=1000000;
    long long atomParsed=0;
    while (1) {
        Snapshot *s=new Snapshot(fXyz);
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
    FILE *fEffectiveRadius=fopen("effectiveRadius.txt","w");
    const int nFrame=trajectory.size();
    const int iBegin=(int)(FRAME_BEGIN*nFrame);
    const int iEnd=nFrame-1;
    fprintf(gnuplot,"set terminal svg;\n");
    fprintf(gnuplot,"set title \'SYS=%s  h=%f A=%f n=%d\';\n",TITLE,Z_HALFDEPTH*2,ANNULUS_AREA,iEnd-iBegin+1);
    vector< pair<double,double> > effectiveRadiusList;
    for (double zCenter=Z_LOW;zCenter<Z_HIGH;zCenter+=Z_INTERVAL) {
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
    fprintf(gnuplot,"set size ratio -1;\n");
    fprintf(gnuplot,"set output \'effectiveRadius.svg\';\n");
    fprintf(gnuplot,"plot \'effectiveRadius.txt\' u 2:1 t \'\',");
    {
        pair<double,double> circleParam=fitCircle(effectiveRadiusList);
        double h=circleParam.first;
        double r=circleParam.second;
        fprintf(gnuplot,"sqrt(%f*%f-x*x)%s%f t \'cos(theta)=%f\';\n",r,r,h>=0?"+":"",h,-h/r);
    }
    for (int i=0;i<trajectory.size();i++) {
        delete trajectory[i];
    }
    return 0;
}

