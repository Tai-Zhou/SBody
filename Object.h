#ifndef _OBJECT_H
#define _OBJECT_H

namespace SBody {
	namespace Object {
		class object {
		  public:
			virtual int hit(const double current[], const double last[]);
		};
		class star : public object {
		  private:
			const int fixed;
			const double radius;
			double pos[8];

		  public:
			star(double radius, const double position[], int fixed = 0);
			int hit(const double current[], const double last[]);
		};
		class disk : public object {
		  private:
			const double innerRadius;
			const double outerRadius;
			const double halfAngle;

		  public:
			disk(double innerRadius, double outerRadius, double halfAngle);
			virtual int hit(const double current[], const double last[]);
		};
		class torus : public object {
		  private:
			const double majorRadius;
			const double minorRadius;

		  public:
			torus(double majorRadius, double minorRadius);
			virtual int hit(const double current[], const double last[]);
		};
	} // namespace Object
} // namespace SBody

#endif
