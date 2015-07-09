//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Implementation of the krovak (Krovak) projection.
//			Definition: http://www.ihsenergy.com/epsg/guid7.html#1.4.3
// Author:	Thomas Flemming, tf@ttqv.com
//
//*****************************************************************************
// Copyright (c) 2001, Thomas Flemming, tf@ttqv.com
// Copyright (c) 2008-2011 by the Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//*****************************************************************************

using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_krovak : PJ
	{
		protected double C_x;

		public override string Name { get { return "krovak"; } }
		public override string DescriptionName { get { return "Krovak"; } }
		public override string DescriptionType { get { return "PCyl, Ell"; } }
		public override string DescriptionParameters { get { return "lat_ts="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				if(C_x!=0) ret.AppendFormat(nc, " +lat_ts={0}", C_x*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		// NOTES:	According to EPSG the full Krovak projection method should have
		//			the following parameters. Within PROJ.4 the azimuth, and pseudo
		//			standard parallel are hardcoded in the algorithm and can't be
		//			altered from outside. The others all have defaults to match the
		//			common usage with Krovak projection.
		//
		// lat_0 = latitude of centre of the projection
		//
		// lon_0 = longitude of centre of the projection
		//
		// ** = azimuth (true) of the centre line passing through the centre of the projection
		//
		// ** = latitude of pseudo standard parallel
		//
		// k = scale factor on the pseudo standard parallel
		//
		// x_0 = False Easting of the centre of the projection at the apex of the cone
		//
		// y_0 = False Northing of the centre of the projection at the apex of the cone

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			// calculate xy from lat/lon

			// Constants, identical to inverse transform function
			double s45=0.785398163397448;	// 45deg
			double s90=2*s45;
			double fi0=phi0;				// Latitude of projection centre 49deg 30'

			// Ellipsoid is used as Parameter in pj_fwd.cs and pj_inv.cs, therefore a must
			// be set to 1 here.
			// Ellipsoid Bessel 1841 a = 6377397.155m 1/f = 299.1528128,
			// e2=0.006674372230614;
			double a=1; //6377397.155;
			//e2=P.es;
			double e2=0.006674372230614;
			double e=Math.Sqrt(e2);

			double alfa=Math.Sqrt(1.0+(e2*Math.Pow(Math.Cos(fi0), 4))/(1.0-e2));
			double uq=1.04216856380474;		// DU(2, 59, 42, 42.69689)
			double u0=Math.Asin(Math.Sin(fi0)/alfa);
			double g=Math.Pow((1.0+e*Math.Sin(fi0))/(1.0-e*Math.Sin(fi0)), alfa*e/2.0);
			double k=Math.Tan(u0/2.0+s45)/Math.Pow(Math.Tan(fi0/2.0+s45), alfa)*g;
			double k1=k0;
			double n0=a*Math.Sqrt(1.0-e2)/(1.0-e2*Math.Pow(Math.Sin(fi0), 2));
			double s0=1.37008346281555;		// Latitude of pseudo standard parallel 78deg 30'00" N
			double n=Math.Sin(s0);
			double ro0=k1*n0/Math.Tan(s0);
			double ad=s90-uq;

			// Transformation
			double gfi=Math.Pow(((1.0+e*Math.Sin(lp.phi))/(1.0-e*Math.Sin(lp.phi))), (alfa*e/2.0));
			double u=2.0*(Math.Atan(k*Math.Pow(Math.Tan(lp.phi/2.0+s45), alfa)/gfi)-s45);
			double deltav=-lp.lam*alfa;
			double s=Math.Asin(Math.Cos(ad)*Math.Sin(u)+Math.Sin(ad)*Math.Cos(u)*Math.Cos(deltav));
			double d=Math.Asin(Math.Cos(u)*Math.Sin(deltav)/Math.Cos(s));
			double eps=n*d;
			double ro=ro0*Math.Pow(Math.Tan(s0/2.0+s45), n)/Math.Pow(Math.Tan(s/2.0+s45), n);

			// x and y are reverted!
			xy.y=ro*Math.Cos(eps)/a;
			xy.x=ro*Math.Sin(eps)/a;

			if(!Proj.pj_param_t(ctx, parameters, "czech"))
			{
				xy.y*=-1.0;
				xy.x*=-1.0;
			}

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			// calculate lat/lon from xy

			// Constants, identisch wie in der Umkehrfunktion
			double s45=0.785398163397448;		// 45deg
			double s90=2*s45;
			double fi0=phi0;					// Latitude of projection centre 49deg 30'

			// Ellipsoid is used as Parameter in for.c and inv.c, therefore a must 
			// be set to 1 here.
			// Ellipsoid Bessel 1841 a = 6377397.155m 1/f = 299.1528128,
			double a=1; // 6377397.155;
			//e2=P.es;
			double e2=0.006674372230614;
			double e=Math.Sqrt(e2);

			double alfa=Math.Sqrt(1.0+(e2*Math.Pow(Math.Cos(fi0), 4))/(1.0-e2));
			double uq=1.04216856380474;		// DU(2, 59, 42, 42.69689)
			double u0=Math.Asin(Math.Sin(fi0)/alfa);
			double g=Math.Pow((1.0+e*Math.Sin(fi0))/(1.0-e*Math.Sin(fi0)), alfa*e/2.0);
			double k=Math.Tan(u0/2.0+s45)/Math.Pow(Math.Tan(fi0/2.0+s45), alfa)*g;
			double k1=k0;
			double n0=a*Math.Sqrt(1.0-e2)/(1.0-e2*Math.Pow(Math.Sin(fi0), 2));
			double s0=1.37008346281555;		// Latitude of pseudo standard parallel 78deg 30'00" N
			double n=Math.Sin(s0);
			double ro0=k1*n0/Math.Tan(s0);
			double ad=s90-uq;

			// Transformation
			// revert y, x
			double xy0=xy.x;
			xy.x=xy.y;
			xy.y=xy0;

			if(!Proj.pj_param_t(ctx, parameters, "czech"))
			{
				xy.x*=-1.0;
				xy.y*=-1.0;
			}

			double ro=Math.Sqrt(xy.x*xy.x+xy.y*xy.y);
			double eps=Math.Atan2(xy.y, xy.x);
			double d=eps/Math.Sin(s0);
			double s=2.0*(Math.Atan(Math.Pow(ro0/ro, 1.0/n)*Math.Tan(s0/2.0+s45))-s45);

			double u=Math.Asin(Math.Cos(ad)*Math.Sin(s)-Math.Sin(ad)*Math.Cos(s)*Math.Cos(d));
			double deltav=Math.Asin(Math.Cos(s)*Math.Sin(d)/Math.Cos(u));

			lp.lam=lam0-deltav/alfa;

			// ITERATION FOR lp.phi
			double fi1=u;
			bool ok=false;
			do
			{
				lp.phi=2.0*(Math.Atan(Math.Pow(k, -1.0/alfa)*Math.Pow(Math.Tan(u/2.0+s45), 1.0/alfa)*
					Math.Pow((1.0+e*Math.Sin(fi1))/(1.0-e*Math.Sin(fi1)), e/2.0))-s45);
				if(Math.Abs(fi1-lp.phi)<0.000000000000001) ok=true;
				fi1=lp.phi;

			} while(!ok);

			lp.lam-=lam0;

			return lp;
		}

		public override PJ Init()
		{
			// read some Parameters,
			// here Latitude Truescale

			double ts=Proj.pj_param_r(ctx, parameters, "lat_ts");
			C_x=ts;

			// we want Bessel as fixed ellipsoid
			a=6377397.155;
			es=0.006674372230614;
			e=Math.Sqrt(es);

			// if latitude of projection center is not set, use 49d30'N
			if(!Proj.pj_param_t(ctx, parameters, "lat_0")) phi0=0.863937979737193;

			// if center long is not set use 42d30'E of Ferro - 17d40' for Ferro
			// that will correspond to using longitudes relative to greenwich
			// as input and output, instead of lat/long relative to Ferro
			if(!Proj.pj_param_t(ctx, parameters, "lon_0")) lam0=0.7417649320975901-0.308341501185665;

			// if scale not set default to 0.9999
			if(!Proj.pj_param_t(ctx, parameters, "k")) k0=0.9999;

			// always the same
			inv=e_inverse;
			fwd=e_forward;

			return this;
		}
	}
}
