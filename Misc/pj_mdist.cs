// libproj -- library of cartographic projections
//
// Copyright (c) 2003, 2006 Gerald I. Evenden
// Copyright (c) 2009-2011 by the Authors
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// Computes distance from equator along the meridian to latitude phi
// and inverse on unit ellipsoid.
// Precision commensurate with double precision.

using System;

namespace Free.Ports.Proj4
{
	class MDIST
	{
		const int MAX_ITER_mdist=20;
		const double TOL14=1.0e-14;

		public int nb;
		public double es;
		public double E;
		public double[] b;

		public MDIST(double es)
		{
			double[] E=new double[MAX_ITER_mdist];

			// generate E(e^2) and its terms E[]
			double ens=es;
			double numf=1.0;
			double twon1=1.0;
			double denfi=1.0;
			double denf=1.0;
			double twon=4.0;
			double Es=1.0;
			double El=1.0;
			E[0]=1.0;

			int i=1;
			for(; i<MAX_ITER_mdist; ++i)
			{
				numf*=(twon1*twon1);
				double den=twon*denf*denf*twon1;
				double T=numf/den;
				E[i]=T*ens;
				Es-=E[i];
				ens*=es;
				twon*=4.0;
				denfi++;
				denf*=denfi;
				twon1+=2.0;
				if(Es==El) break; // jump out if no change
				El=Es;
			}

			this.b=new double[i+1];

			this.nb=i-1;
			this.es=es;
			this.E=Es;

			// generate b_n coefficients--note: collapse with prefix ratios
			this.b[0]=Es=1.0-Es;
			numf=denf=1.0;
			double numfi=2.0;
			denfi=3.0;
			for(int j=1; j<i; ++j)
			{
				Es-=E[j];
				numf*=numfi;
				denf*=denfi;
				this.b[j]=Es*numf/denf;
				numfi+=2.0;
				denfi+=2.0;
			}
		}

		public double proj_mdist(double phi, double sphi, double cphi)
		{
			double sc=sphi*cphi;
			double sphi2=sphi*sphi;
			double D=phi*E-es*sc/Math.Sqrt(1.0-es*sphi2);
			int i=nb;
			double sum=b[i];
			while(i-->0) sum=b[i]+sphi2*sum;

			return D+sc*sum;
		}

		public double proj_inv_mdist(projCtx ctx, double dist)
		{
			double k=1.0/(1.0-es);
			int i=MAX_ITER_mdist;
			double phi=dist;
			while(i-->0)
			{
				double s=Math.Sin(phi);
				double t=1.0-es*s*s;
				phi-=t=(proj_mdist(phi, s, Math.Cos(phi))-dist)*(t*Math.Sqrt(t))*k;
				if(Math.Abs(t)<TOL14) return phi; // that is no change
			}

			// convergence failed
			Proj.pj_ctx_set_errno(ctx, -17);
			return phi;
		}
	}
}
