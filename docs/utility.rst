.. index::
   Utility

=======
Utility
=======

.. doxygenvariable:: SBody::absolute_accuracy
   :no-link:
.. doxygenvariable:: SBody::relative_accuracy
   :no-link:
.. doxygenvariable:: SBody::sample_number
   :no-link:
.. doxygenvariable:: SBody::epsilon
   :no-link:
.. doxygenvariable:: SBody::sin_epsilon
   :no-link:
.. doxygenvariable:: SBody::cos_epsilon
   :no-link:
.. doxygenvariable:: SBody::M_2PI
   :no-link:
.. doxygenvariable:: SBody::M_PI2
   :no-link:
.. doxygenvariable:: SBody::M_SQRT27
   :no-link:

Integrator
==========

.. doxygenclass:: SBody::Integrator
   :members:
   :no-link:

Solver
======

.. doxygenclass:: SBody::Solver
   :members:
   :no-link:

FunctionSolver
--------------

.. doxygenclass:: SBody::FunctionSolver
   :members:
   :no-link:

DerivativeSolver
----------------

.. doxygenclass:: SBody::DerivativeSolver
   :members:
   :no-link:


MultiSolver
===========

.. doxygenclass:: SBody::MultiSolver
   :members:
   :no-link:

MultiFunctionSolver
-------------------

.. doxygenclass:: SBody::MultiFunctionSolver
   :members:
   :no-link:

MultiDerivativeSolver
---------------------

.. doxygenclass:: SBody::MultiDerivativeSolver
   :members:
   :no-link:

GslBlock
========

.. doxygenclass:: SBody::GslBlock
   :members:
   :no-link:

Math
====

.. doxygenfunction:: SBody::Dot(const double[], const double[], size_t);
   :no-link:

.. doxygenfunction:: SBody::Dot(const double[], size_t);
   :no-link:

.. doxygenfunction:: SBody::Norm
   :no-link:

.. doxygenfunction:: SBody::Cross
   :no-link:

.. doxygenfunction:: SBody::DotCross
   :no-link:

.. doxygenfunction:: SBody::TriangleArea(double, double, double)
   :no-link:

.. doxygenfunction:: SBody::TriangleArea(const double[], const double[], const double[])
   :no-link:

.. doxygenfunction:: SBody::RotateAroundAxis
   :no-link:

.. doxygenfunction:: SBody::CartesianToSpherical(const double[], double[], size_t)
   :no-link:

.. doxygenfunction:: SBody::CartesianToSpherical(double[], size_t)
   :no-link:

.. doxygenfunction:: SBody::SphericalToCartesian(const double[], double[], size_t)
   :no-link:

.. doxygenfunction:: SBody::SphericalToCartesian(double[], size_t)
   :no-link:

.. doxygenfunction:: SBody::OppositeSign
   :no-link:

.. doxygenfunction:: SBody::MapTheta
   :no-link:

.. doxygenfunction:: SBody::ModBy2Pi
   :no-link:

.. doxygenfunction:: SBody::PhiDifference
   :no-link:

.. doxygenfunction:: SBody::LinearInterpolation(double, double, double, double, double)
   :no-link:

.. doxygenfunction:: SBody::LinearInterpolation(double, double, double, const double[], const double[], double[], size_t)
   :no-link:

.. doxygenfunction:: SBody::InterpolateSphericalPositionToCartesian
   :no-link:

Elliptic Integrals
------------------

Calculate the elliptic integrals using the Legendre forms. Further information can be found in `GSL Manual <https://www.gnu.org/software/gsl/doc/html/specfunc.html#elliptic-integrals>`_, `Carlson (1988) <https://www.ams.org/mcom/1988-51-183/S0025-5718-1988-0942154-7/>`_, `Carlson (1989) <https://www.ams.org/mcom/1989-53-187/S0025-5718-1989-0969482-4/>`_, `Carlson (1991) <https://www.ams.org/mcom/1991-56-193/S0025-5718-1991-1052087-6/>`_, and `Carlson (1992) <https://www.ams.org/mcom/1992-59-199/S0025-5718-1992-1134720-4/>`_.

.. doxygenfunction:: SBody::EllipticIntegral
   :no-link:

.. doxygenfunction:: SBody::EllipticIntegral2Complex
   :no-link:

.. doxygenfunction:: SBody::EllipticIntegral4Complex
   :no-link:

.. doxygenfunction:: SBody::Carlson_RC
   :no-link:

.. doxygenfunction:: SBody::Carlson_RJ
   :no-link:
