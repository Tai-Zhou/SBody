#ifndef _OBJECT_H
#define _OBJECT_H

namespace SBody {
	namespace Object {
		class object {
		  public:
			virtual int hit(const double current[], const double last[]);
		};
		class star : public object {
		  protected:
			const int fixed;
			const double radius;
			double pos[8];

		  public:
			star(double radius, const double position[], int fixed = 0);
			int hit(const double current[], const double last[]);
		};
		class disk : public object {
		  protected:
			const double innerRadius;
			const double outerRadius;

		  public:
			disk(double innerRadius, double outerRadius);
			int hit(const double current[], const double last[]);
		};
		class thickDisk : public disk {
		  protected:
			const double halfAngle;

		  public:
			thickDisk(double innerRadius, double outerRadius, double halfAngle);
			int hit(const double current[], const double last[]);
		};
		class torus : public object {
		  protected:
			const double majorRadius;
			const double minorRadius;

		  public:
			torus(double majorRadius, double minorRadius);
			int hit(const double current[], const double last[]);
		};
	} // namespace Object
} // namespace SBody

#endif
