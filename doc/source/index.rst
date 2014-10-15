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
selection or to calculate photometric redshifts for quasars, using
any combination of SDSS optical, GALEX, ultraviolet, UKIDSS near-IR,
and WISE mid-IR photometry

Contents:

	:ref:`Introduction <intro>`

	:ref:`XDQSOZ photometric quasar catalog <cat>`

	:ref:`IDL code <idl>`

	:ref:`Acknowledging XDQSO <ack>`



.. _intro:

Introduction
-------------

To download the code use either

.. code-block:: none

   wget -qO- https://github.com/xdqso/xdqso/archive/v0.6.tar.gz | tar xvz

or

.. code-block:: none

   git clone https://github.com/xdqso/xdqso.git


Installation only requires you to set the environment variable
``XDQSODATA`` to the ``data`` directory of the distribution. EvilUPS
setup is available.

Code is available either as ``xdqso`` or as ``xdqsoz``. For most
purposes you will want to use the ``xdqsoz`` functions: these allow
you to calculate photometric redshifts and quasar probabilities for
arbitrary redshift ranges. The ``xdqsoz`` routines are the functions
used to create the photometric quasar catalog. If you want to mimick
*SDSS-III's BOSS* quasar target selection or easily calculate
probabilites in broad redshift bins you want to use the
``xdqso`` functions instead.

The functions ``xdqso_calculate_prob`` and ``xdqsoz_calculate_prob``
calculate photometric quasar probabilities. The former can only do
this in three redshift ranges (*z* :math:`<` 2.2; 2.2 :math:`\leq` *z*
:math:`\leq` 3.5; z :math:`>` 3.5), while the latter accomodates arbitrary
redshift ranges.

Photometric redshifts can be calculated using the ``xdqsoz_zpdf`` and
``xdqsoz_eval_zpdf`` functions. The former prepares the parameters of
the one-dimensional redshift PDF for individual objects, the latter
then allows you to evaluate this PDF.  The function ``xdqsoz_calculate_prob_andz``
wraps this functionality into the quasar probability calculation, including
the redshift PDF sampled in increments of 0.01 in the output structure.

One can also use xdqsoz_marginalize_colorzprob to integrate the redshift
PDF over arbitrary redshift ranges:

.. code-block:: none

   out= xdqsoz_marginalize_colorzprob(zmin,zmax,flux,flux_ivar,norm=totlike)

Input is dereddened psfflux and psfflux_ivar (to deredden you can use
the functions xdqso_sdss_deredden and xdqso_sdss_deredden_error) and a
min and max redshift; output is the marginalized likelihood
(marginalized over redshift). Setting norm=totlike returns the total
quasar likelihood. If you then multiply the 'pqso' from the
photometric quasar catalog below by the ratio of out and totlike you
get the desired redshift-range probability (since the prior and the
denominator do not change). When calculating quasar probabilities in
many bins this is much faster than calling ``xdqsoz_calculate_prob``
repeatedly because you do not recalculate star likelihoods and priors
each time.

.. _cat:

XDQSOz photometric quasar catalog
---------------------------------

The original version of the SDSS-XDQSO DR8 photometric quasar
catalog that does not include WISE information is available at

http://cosmo.nyu.edu/~jb2777/qsocat/xdqsoz_pqso0.5_imag21.5-nobadu.fits.gz

http://cosmo.nyu.edu/~jb2777/qsocat/README_pqso0.5_imag21.5-nobadu

This catalog is based on the same principle as the XDQSO method for
BOSS quasar selection, but uses a slightly different algorithm
(*XDQSOz*) for calculating quasar probabilities that also permits us
to obtain photometric redshifts; it also allows quasar probabilities
to be calculated quickly for arbitrary redshift ranges (see the
accompanying code below).

The original catalog is a simple cut on P(quasar) > 0.5 for all
objects that pass the BOSS quasar selection flag cuts, limited further
to i_0 < 21.5 mag and with some bad u-columns in the SDSS imaging data
masked. We have performed some first tests of the clustering of the
objects in the catalog, which shows that the level of stellar
contamination is small (< 10%), but we have yet to break this up by
redshift range, etc., and perform further tests. So exercise caution
when using the catalog (especially at low Galactic latitude, since the
SEGUE stripes are included), and please let us know if you find any
problems.

An updated version of the catalog is available at 

[link to be inserted - please contact code owners for temporary link]

which includes updated probabilities incorporporating WISE fluxes, and
photometric redshift PDFs for all objects with P(quasar) > 0.2.  Like the
first catalog, it only includes objects that pass the BOSS quasar selection 
flag cuts and objects with i_0 < 21.5.  There are tags to indicate if an object
falls within the SDSS bright star mask, a region of bad SDSS photometry, 
an area with bad u-columns, or near contaminated WISE data.  The same precautions
as above apply to the new catalog.


.. _idl:


IDL code
--------

Contents:

	:ref:`xdqso_calculate_prob <idl_xdqso_calculate_prob>`

	:ref:`xdqsoz_calculate_prob <idl_xdqsoz_calculate_prob>`

	:ref:`xdqsoz_eval_zpdf <idl_xdqsoz_eval_zpdf>`

	:ref:`xdqsoz_marginalize_colorzprob <idl_xdqsoz_marginalize_colorzprob>`

	:ref:`xdqsoz_peaks <idl_xdqsoz_peaks>`

	:ref:`xdqsoz_qso_track <idl_xdqsoz_qso_track>`

	:ref:`xdqsoz_zpdf <idl_xdqsoz_zpdf>`

	:ref:`xdqsoz_calculate_prob_andz <idl_xdqsoz_calculate_prob_andz>`

.. _idl_xdqso_calculate_prob:

**xdqso_calculate_prob** (in,/dereddened)

	*calculate the extreme-deconvolution XDQSO QSO probability*

	Input:

		in - structure containing PSFFLUX, PSFFLUX_IVAR, EXTINCTION

	Keywords:

		dereddened - psfflux, and psfflux_ivar is already dereddened

		galex - GALEX fluxes are included in input structure, with tags NUV, FUV, 
		           NUV_ivar, and FUV_ivar.  GALEX fluxes are in nanomaggies

		ukidss - UKIDSS fluxes are included in input structure, with tags APERCSIFLUX3_Y,
                            APERCSIFLUX3_J,  APERCSIFLUX3_H,  APERCSIFLUX3_K, APERCSIFLUX3ERR_Y, 
			    APERCSIFLUX3ERR_J, APERCSIFLUX3ERR_H, APERCSIFLUX3ERR_K.  Fluxes/errors are in SI units.

		wise - WISE fluxes are included in input structure, with tags w1_nanomaggies,
		          w2_nanomaggies, w1_nanomaggies_ivar, w2_nanomaggies_ivar.  Fluxes are in
		          Vega nanomaggies.

	Output:

		structure containing pqso, ... (see XDQSO catalog description)
			 

	History:

		010-04-30 - Written - Bovy (NYU)

		2014-04-02 - Added WISE, GALEX, UKIDSS - DiPompeo (UWyo)


.. _idl_xdqsoz_calculate_prob:

**xdqsoz_calculate_prob** (in,zmin,zmax,/dereddened,/galex,/ukidss)

	*calculate the extreme-deconvolution probability ratio, marginalizing over an arbitrary redshift range*

	Input:

		in - structure containing PSFFLUX, PSFFLUX_IVAR, EXTINCTION

		zmin, zmax - lower, upper bound of redshift interval

	Keywords:

		dereddened  - psfflux, and psfflux_ivar are already dereddened

		galex - GALEX fluxes are included in input structure, with tags NUV, FUV, 
		          NUV_ivar, and FUV_ivar.  GALEX fluxes are in nanomaggies

		ukidss - UKIDSS fluxes are included in input structure, with tags APERCSIFLUX3_Y,
                  		APERCSIFLUX3_J,  APERCSIFLUX3_H,  APERCSIFLUX3_K, APERCSIFLUX3ERR_Y,   
          	                APERCSIFLUX3ERR_J, APERCSIFLUX3ERR_H, APERCSIFLUX3ERR_K.  Fluxes/errors are in SI units.

		wise - WISE fluxes are included in input structure, with tags w1_nanomaggies,
		           w2_nanomaggies, w1_nanomaggies_ivar, w2_nanomaggies_ivar.  Fluxes are in
		           Vega nanomaggies.


	Output:

		out - structure containing pqso, ...

	History:

		2010-04-30 - Written - Bovy (NYU)

		2010-05-29 - Added Galex - Bovy

		2010-10-30 - Added UKIDSS - Bovy

		2014-03-31 - Added WISE - DiPompeo (UWyo)


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


.. _idl_xdqsoz_marginalize_colorzprob:

**xdqsoz_marginalize_colorzprob** (zmin,zmax,flux,flux_ivar,/galex,/ukidss,norm=norm,/log)

	*marginalize the probability of a relative flux + redshift (not a color) over redshift*

	Input:

		zmin, zmax - redshift

		flux - [nfluxes] or [nfluxes,ndata] array of fluxes

		flux_ivar - [nfluxes] or [nfluxes,ndata] array of flux_ivars
	
	Keywords:

		galex - use GALEX fits

		ukidss - use UKIDSS fits
		
		wise - use WISE fits

		log - calculate log

	Output:

		number or array of probabilities

	Optional Output:
	
		norm - normalization factor (likelihood marginalized over redshift 0 to infinity)

	History:

		2011-01-16 - Written - Bovy (NYU)

		2014-03-31 - Added WISE - DiPompeo (UWyo)


.. _idl_xdqsoz_peaks:

**xdqsoz_peaks** (flux,flux_ivar,nzs=nzs,peak_threshold=peak_threshold,/galex,/ukidss,/plot,peakz=peakz,xdqsoz=xdqsoz)

        *calculate the number of peaks of a zpdf as well as the MAP z*

	Input:

		flux - dereddened flux

		flux_ivar - dereddened flux_ivar

	Optional Input:

		 nzs - number of points to sample the PDF at

		 peak_threshold - threshold for defining a peak (contiguous region with p above peak_threshold)

	Keywords:

		galex - use GALEX fits

		ukidss - use UKIDSS fits
		
		wise - use WISE fits

		plot - make QS plot

	Output:
	
		number of peaks
	
	Optional Output:

		 peakz - MAP z

		 xdqsoz - structure containing {peakz,peakprob,peakfwhm,otherz,otherprob,otherfwhm} for all peaks

	History:

		2011-01-18 - Written - Bovy (NYU)

		2014-03-31 - Added WISE - DiPompeo (UWyo)


.. _idl_xdqsoz_qso_track:

**xdqsoz_qso_track** (z,i=i,/galex,/ukidss)

        *calculate the mean quasar locus*

	Input:

		z - redshift or array of redshifts [N]

	Optional Input:

	      i= dereddened i-band magnitude

	Keywords: 

		galex - use GALEX fits

		ukidss - use UKIDSS fits

		wise - use WISE fits

	Output:

		mags[ndim,N] - array of apparent magnitudes (ugriz[NUV,FUV,YJHK])

	History:

		2011-04-01 - Written - Bovy (NYU)

		2014-04-02 - Added WISE - DiPompeo (UWyo)


.. _idl_xdqsoz_zpdf:

**xdqsoz_zpdf**, flux, flux_ivar, /galex, /ukidss, zmean=zmean, zcovar=zcovar, zamp=zamp

	*calculate the photometric redshift pdf using XDQSOz*

	Input:

		flux - [nfluxes] or [nfluxes,ndata] array of fluxes
		
		flux_ivar - [nfluxes] or [nfluxes,ndata] array of flux_ivars

	Keywords:

		galex - use GALEX fits
		
		ukidss - use UKIDSS fits

		wise - use WISE fits
	
	Output:

		zmean - [ngauss,ndata] array of means
		
		zcovar - [ngauss,ndata] array of covars
		
		zamp - [ngauss,ndata] array of amplitudes

	History:
	  
		2011-01-18 - Written - Bovy (NYU)

		2014-04-02 - Added WISE - DiPompeo (UWyo)


.. _idl_xdqsoz_calculate_prob_andz:

**xdqsoz_calculate_prob** (in,zmin,zmax,/dereddened,/galex,/ukidss)

	*The same as xdqsoz_calculate_prob, with xdqsoz_zpdf wrapped in to simultaneously calculate z PDF*

	Input:

		in - structure containing PSFFLUX, PSFFLUX_IVAR, EXTINCTION

		zmin, zmax - lower, upper bound of redshift interval

	Keywords:

		dereddened  - psfflux, and psfflux_ivar are already dereddened

		galex - GALEX fluxes are included in psfflux, psfflux_ivar, and extinction; use them

		ukidss - use UKIDSS (like /galex)

		wise - use WISE (like /galex)

	Output:

		out - structure containing pqso, ... , z array from zmin to zmax in 0.01 increments,
		      z PDF at each value of z.

	History:

		2014-03-31 - Written - DiPompeo (UWyo)





.. _ack:

Acknowledging XDQSO
--------------------

Please cite the relevant papers among the following:

       BOSS CORE target selection paper (also cite `Ross et al. 2011 <http://adsabs.harvard.edu/abs/2011arXiv1105.0606R>`_): *Think outside the color box: probabilistic target selection and the SDSS-XDQSO quasar targeting catalog*, Bovy, J., et al., 2010, ApJ, **729**, 141 `[ApJ] <http://dx.doi.org/10.1088/0004-637X/729/2/141>`_ `[ADS] <http://adsabs.harvard.edu/abs/2011ApJ....729..141B>`_


       Photometric redshifts: *Photometric redshifts and quasar probabilities from a single, data-driven generative model*, Bovy, J., et al., 2011, ApJ, **749**, 41 `[ApJ] <http://dx.doi.org/10.1088/0004-637X/749/1/41>`_ `[ADS] <http://adsabs.harvard.edu/abs/2012ApJ...749...41B>`_


       *Incorporating WISE Photometry into Quasar Probabilities and Photometric Redshift Estimation With XDQSOz*, DiPompeo,M.A., et al., 2014, in preparation


       Catalog paper: *The SDSS-XDQSO photometric quasar catalog*, Myers, A. D., et al., 2015, in preparation


       XD methodology paper: *Extreme deconvolution: inferring complete distribution functions from noisy, heterogeneous and incomplete observations*, Bovy, J., Hogg, D. W., & Roweis, S. T., 2011, AOAS, **5**, 2B, 1657 `[AOAS] <http://dx.doi.org/10.1214/10-AOAS439>`_ `[ADS] <http://adsabs.harvard.edu/abs/2011AnApS...5.1657B>`_


..
	Indices and tables
	==================

	* :ref:`genindex`
	* :ref:`modindex`
	* :ref:`search`

