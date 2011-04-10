#ifndef _IRTKBEP_H

#define _IRTKBEP_H

#include <irtkRegistration.h>
#include <irtkGaussianBlurring.h>

typedef struct {
int cx,cy;
int x1,y1;
int x2,y2;
double dx,dy;
double sx,sy;
}Sector;


class irtkBep : public irtkObject
{
private:
	irtkGreyImage target,source;
	irtkPointSet landmarks;
	irtkPointSet olandmarks;
	irtkPointSet centerpoints;
	int tip,aa,am,mb,bottom;
	char *outputfilename;
	char *maxfilename;
	char *osequence;
	char *lsequence;
	irtkRealImage threshold,region;
	irtkRealImage wallthreshold;
	irtkRealImage segmentation,compare;
	double wallaverage;

public:

	irtkBep();

	virtual void SetCompare(irtkRealImage& compare);

	virtual void SetSegmentation(irtkRealImage& segmentation);

	virtual void SetInput(irtkGreyImage& target, irtkGreyImage& source, irtkRealImage& threshold, char* osequence, char* lsequence);

	virtual void SetInput(irtkGreyImage& target, irtkRealImage& threshold, char* osequence);

	virtual void SetLandmarks(irtkPointSet& olandmarks, int tip, int aa, int am, int mb, int bottom);

	virtual void SetOutput(char *outputfilename, char *maxfilename);

	virtual void Initialize();

	virtual void Finalize();

    virtual void Bullseyeplot(int average = 1);

	virtual void GenerateSegmentation(char* surfacefilename);
	
	virtual void Compareoverall();

	virtual void Comparedifference();

	virtual void Compareradius();

	virtual void AnalysisRadius(int prefix = 0);

	virtual void AnalysisMotion(int prefix = 3);

	virtual void AnalysisStrain(int prefix = 3, int mode = 0);

	virtual void EvaluateWallaverage();

	virtual void EvaluateWallThreshold();

	double **bep;

private:

	void evaluatecenter(Sector& s, int k){
		s.cx = s.sx + k*s.dx;
		s.cy = s.sy + k*s.dy;
	}

	double area(Sector& abc){
		double area = (abc.cx-abc.x2)*(abc.y1-abc.y2) - (abc.cy-abc.y2)*(abc.x1-abc.x2);
		return abs(area);
	}

	bool isinside(int i, int j, Sector& abc){
		Sector pbc,apc,abp;
		pbc = abc;
		apc = abc;
		abp = abc;
		pbc.cx = i; pbc.cy = j;
		apc.x1 = i; apc.y1 = j;
		abp.x2 = i; abp.y2 = j;
		if (abs(area(abc) - area(pbc) - area(apc)
			- area(abp)) < 0.1){
				return 1;
		}else
			return 0;	
	}

	double stintegral(irtkRealImage& inte, Sector& sector, int minz, int maxz, int number, int t, int average = 1){
		int i,j,k;
		double integral = 0;
		double count = 0;
		double k1,k2,kc;
		double r = 0,rmin = 9999, rmax = 0, rtmp = 0;
		if(average != -1 && average != 1){
			for ( k = minz; k < maxz; k++){
				for ( j = 0; j< inte.GetY(); j++){
					for ( i = 0; i<inte.GetX(); i++){
						if(round(segmentation.GetAsDouble(i,j,k)) == number+1){
							integral += inte.GetAsDouble(i,j,k,t);
							rtmp = sqrt(pow((i-sector.cx)*inte.GetXSize(),2)
								+pow((j-sector.cy)*inte.GetYSize(),2));
							r += rtmp;
							if(rtmp > rmax) rmax = rtmp;
							if(rtmp < rmin) rmin = rtmp;
							count = count + 1;
						}
					}
				}
			}
			r = r/count;
			if(rmax - r > r - rmin)
				rmax = 2*r - rmin;
			else
				rmin = 2*r - rmax;
			if(average == 0){	
				// i = size*n/a*s
				k1 = atan2((double)sector.cy - sector.y1,(double)sector.cx - sector.x1);
				k2 = atan2((double)sector.cy - sector.y2,(double)sector.cx - sector.x2);
				kc = abs(k2-k1);
				if(kc > M_PI)
					kc = 2*M_PI - kc;
				integral = integral*inte.GetXSize()*inte.GetYSize()
					/(kc*(maxz-minz));
				integral = integral/r;
			}
		}else if(average == 1){
			for ( k = 0; k < inte.GetZ(); k++){
				for ( j = 0; j< inte.GetY(); j++){
					for ( i = 0; i<inte.GetX(); i++){
						if(round(segmentation.GetAsDouble(i,j,k)) == number+1
							&& inte.GetAsDouble(i,j,k,t) > 0){
							integral += inte.GetAsDouble(i,j,k,t);
							count++;
						}
					}
				}
			}
			integral = integral/count;
		}else{
			for ( k = 0; k < inte.GetZ(); k++){
				for ( j = 0; j< inte.GetY(); j++){
					for ( i = 0; i<inte.GetX(); i++){
						if(round(segmentation.GetAsDouble(i,j,k)) == number+1)
						{
							if(wallthreshold.GetAsDouble(i,j,k) == 2)
								integral++;
							count ++;
						}
					}
				}
			}
			integral = integral*100/count;
		}
		return integral;

	}

	double evaluatesector(irtkRealImage& inte, Sector& sector, irtkPoint &p1,irtkPoint &p2, int minz, int maxz, int number, int t, int average = 1){
		sector.x1 = p1._x; sector.x2 = p2._x;
		sector.y1 = p1._y; sector.y2 = p2._y;
		double bep = stintegral(inte, sector, minz,maxz,number, t, average);
		return bep;
	}

	void segmentsector(irtkRealImage& radius, irtkRealImage& seg, Sector& sector, irtkPoint p1, irtkPoint p2, int min, int max, int label){
		sector.x1 = p1._x;
		sector.y1 = p1._y;
		sector.x2 = p2._x;
		sector.y2 = p2._y;
		double r = 0;
		double count = 0;
		int i,j,k;
		centerpoints(label-1)._x = 0;
		centerpoints(label-1)._y = 0;
		centerpoints(label-1)._z = 0;
		for(k = min; k<max; k++){
			for(j = 0; j<seg.GetY(); j++){
				for(i=0;i<seg.GetX();i++){
					if(isinside(i,j,sector) && threshold.GetAsDouble(i,j,k) == 3){
						seg.PutAsDouble(i,j,k,label);
						r = sqrt(pow((i-sector.cx)*seg.GetXSize(),2)
							+pow((j-sector.cy)*seg.GetYSize(),2));
						radius.PutAsDouble(i,j,k,r);
						centerpoints(label-1)._x += i;
						centerpoints(label-1)._y += j;
						centerpoints(label-1)._z += k;
						count ++;
					}
				}
			}
		}
		centerpoints(label-1)._x = centerpoints(label-1)._x/count;
		centerpoints(label-1)._y = centerpoints(label-1)._y/count;
		centerpoints(label-1)._z = centerpoints(label-1)._z/count;
	}

	void solveequation(double& x1, double& x2, double &y1, double &y2, double Cx, double R, double K, double Cy){
		x1 = Cx - sqrt(pow(R,2)/(pow(K,2) + 1));
		y1 = K*(x1 - Cx) + Cy;

		x2 = Cx + sqrt(pow(R,2)/(pow(K,2) + 1));
		y2 = K*(x2 - Cx) + Cy;
	}

	virtual void analysis(char* sequence,irtkGreyImage& integral);

	virtual void analysis(irtkRealImage& integral);

	virtual void analysis(irtkRealImage& integral, irtkPoint& center);

};

#endif