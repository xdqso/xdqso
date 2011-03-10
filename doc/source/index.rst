.. XDQSO documentation master file, created by
   sphinx-quickstart on Wed Mar  9 16:49:56 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to XDQSO's documentation!
=================================

The XDQSO code accompanies the XDQSO/XDQSOz papers for quasar
classification and photometric redshift estimation. It allows you to
calculate photometric quasar probabilities to mimick `SDSS-III's
BOSS <http://www.sdss3.org/surveys/boss.php>`_ quasar target
selection or to calculate photometric redshifts for quasars.

Contents:

	:ref:`Introduction <intro>`

	:ref:`IDL code <idl>`

	:ref:`Acknowledging XDQSO <ack>`



.. _intro:

Introduction
-------------

To download the code use either

.. code-block:: none

   svn export http://www.sdss3.org/svn/repo/xdqso/trunk xdqso

or

.. code-block:: none

   svn co http://www.sdss3.org/svn/repo/xdqso/trunk xdqso-read-only


Installation only requires you to set the environment variable
``XDQSODATA`` to the ``data`` directory of the distribution (without
the trailing slash). EvilUPS setup is available.


Code is available either as ``xdqso`` or as ``xdqsoz``. For most
purposes you will want to use the ``xdqsoz`` functions: these allow
you to calculate photometric redshifts and quasar probabilities for
arbitrary redshift ranges. The ``xdqsoz`` routines are the functions
used to create the photometric quasar catalog. If you want to mimick
*SDSS-III's BOSS* quasar target selection you want to use the
``xdqso`` functions instead.

The functions ``xdqso_calculate_prob`` and ``xdqsoz_calculate_prob``
calculate photometric quasar probabilities. The former can only do
this in three redshift ranges (*z* :math:`<` 2.2; 2.2 :math:`\leq` *z*
:math:`\leq` 3.5; z :math:`>` 3.5), while the latter accomodates arbitrary
redshift ranges.

Photometric redshifts can be calculate using the ``xdqsoz_zpdf`` and
``xdqsoz_eval_zpdf`` functions. The former prepares the parameters of
the one-dimensional redshift PDF for individual objects, the latter
then allows you to evaluate this PDF.


.. _idl:


IDL code
--------

Contents:

	:ref:`xdqso_calculate_prob <idl_xdqso_calculate_prob>`

	:ref:`xdqsoz_calculate_prob <idl_xdqsoz_calculate_prob>`

	:ref:`xdqsoz_marginalize_colorzprob <idl_xdqsoz_maringalize_colorzprob>`

	:ref:`xdqsoz_eval_zpdf <idl_xdqsoz_eval_zpdf>`

	:ref:`xdqsoz_zpdf <idl_xdqsoz_zpdf>`

.. _idl_xdqso_calculate_prob:

**xdqso_calculate_prob** (in,/dereddened)

	*calculate the extreme-deconvolution XDQSO QSO probability*

	Input:

		in - structure containing PSFFLUX, PSFFLUX_IVAR, EXTINCTION

	Keywords:

		dereddened - psfflux, and psfflux_ivar is already dereddened

	Output:

		structure containing pqso, ... (see XDQSO catalog description)
			 

	History:

		010-04-30 - Written - Bovy (NYU)


.. _idl_xdqsoz_calculate_prob:

**xdqsoz_calculate_prob** (in,zmin,zmax,/dereddened,/galex,/ukidss)

	*calculate the extreme-deconvolution probability ratio, marginalizing over an arbitrary redshift range*

	Input:

		in - structure containing PSFFLUX, PSFFLUX_IVAR, EXTINCTION

		zmin, zmax - lower, upper bound of redshift interval

	Keywords:

		dereddened  - psfflux, and psfflux_ivar are already dereddened

		galex - GALEX fluxes are included in psfflux, psfflux_ivar, and extinction; use them

		ukidss - use UKIDSS (like /galex)

	Output:

		out - structure containing pqso, ...

	History:

		2010-04-30 - Written - Bovy (NYU)


.. _idl_xdqsoz_maringalize_colorzprob:

**xdqsoz_marginalize_colorzprob** (zmin,zmax,flux,flux_ivar,/galex,/ukidss,norm=norm,/log)

	*marginalize the probability of a relative flux + redshift (not a color) over redshift*

	Input:

		zmin, zmax - redshift

		flux - [nfluxes] or [nfluxes,ndata] array of fluxes

		flux_ivar - [nfluxes] or [nfluxes,ndata] array of flux_ivars
	
	Optional Input:
	
		norm - normalization factor (if precomputed by calling this function before with zmin=0.3 and zmax=5.5)

	Keywords:

		galex - use GALEX fits

		ukidss - use UKIDSS

		log - calculate log

	Output:

		number or array of probabilities

	History:

		20111-01-16 - Written - Bovy (NYU)


.. _idl_xdqsoz_eval_zpdf:

**xdqsoz_eval_zpdf** (z,zmean,zcovar,zamp)

	*evaluate the photometric redshift PDF for a given redshift given means, covars, and amps*

	Input:

		z - redshift [nz]
		
		zmean, zcovar, zamp - from :ref:`xdqsoz_zpdf <idl_xdqsoz_zpdf>`

	Output:
	
		p(z)

	History:

		2011-01-18 - Written - Bovy (NYU)


.. _idl_xdqsoz_zpdf:

**xdqsoz_zpdf**, flux, flux_ivar, /galex, /ukidss, zmean=zmean, zcovar=zcovar, zamp=zamp

	*calculate the photometric redshift pdf using XDQSOz*

	Input:

		flux - [nfluxes] or [nfluxes,ndata] array of fluxes
		
		flux_ivar - [nfluxes] or [nfluxes,ndata] array of flux_ivars

	Keywords:

		galex - use GALEX fits
		
		ukidss - use UKIDSS
	
	Output:

		zmean - [ngauss,ndata] array of means
		
		zcovar - [ngauss,ndata] array of covars
		
		zamp - [ngauss,ndata] array of amplitudes

	History:
	  
		2011-01-18 - Written - Bovy (NYU)



.. _ack:

Acknowledging XDQSO
--------------------

Please cite the relevant papers among the following:

       *Think outside the color box: probabilistic target selection and the SDSS-XDQSO quasar targeting catalog*, Bovy, J., et al., 2010, ApJ, **729**, 141 `[ApJ] <http://dx.doi.org/10.1088/0004-637X/729/2/141>`_ `[ADS] <http://adsabs.harvard.edu/abs/2011ApJ....729..141B>`_


       *Photometric redshifts and quasar probabilities over arbitrary redshift ranges*, Bovy, J., et al., 2011, in preparation



       *The SDSS-XDQSO photometric quasar catalog*, Myers, A. D., et al., 2011, in preparation


..
	Indices and tables
	==================

	* :ref:`genindex`
	* :ref:`modindex`
	* :ref:`search`

