//
// libproj -- library of cartographic projections
//
// Copyright (c) 2003, Gerald I. Evenden
// Copyright (c) 2008-2011 by the Authors
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

using System;

namespace Free.Ports.Proj4.Gauss
{
	public class GAUSS
	{
		public double C, K, e, ratexp;

		const double TOL14=1.0e-14;

		static double srat(double esinp, double exp)
		{
			return Math.Pow((1.0-esinp)/(1.0+esinp), exp);
		}

		public static GAUSS pj_gauss_ini(double e, double phi0, out double chi, out double rc)
		{
			chi=rc=0;
			try
			{
				GAUSS en=new GAUSS();

				double es=e*e;
				en.e=e;
				double sphi=Math.Sin(phi0);
				double cphi=Math.Cos(phi0);
				cphi*=cphi;
				rc=Math.Sqrt(1.0-es)/(1.0-es*sphi*sphi);
				en.C=Math.Sqrt(1.0+es*cphi*cphi/(1.0-es));
				chi=Math.Asin(sphi/en.C);
				en.ratexp=0.5*en.C*e;
				en.K=Math.Tan(0.5*chi+Proj.FORTPI)/(Math.Pow(Math.Tan(0.5*phi0+Proj.FORTPI), en.C)*srat(en.e*sphi, en.ratexp));
				return en;
			}
			catch
			{
				return null;
			}
		}

		public static LP pj_gauss(projCtx ctx, LP elp, GAUSS en)
		{
			LP slp;
			slp.phi=2.0*Math.Atan(en.K*Math.Pow(Math.Tan(0.5*elp.phi+Proj.FORTPI), en.C)*srat(en.e*Math.Sin(elp.phi), en.ratexp))-Proj.HALFPI;
			slp.lam=en.C*elp.lam;
			return slp;
		}

		public static LP pj_inv_gauss(projCtx ctx, LP slp, GAUSS en)
		{
			const int MAX_ITER=20;
			LP elp;
			elp.phi=0;

			elp.lam=slp.lam/en.C;
			double num=Math.Pow(Math.Tan(0.5*slp.phi+Proj.FORTPI)/en.K, 1.0/en.C);

			int i=MAX_ITER;
			for(; i>0; i--)
			{
				elp.phi=2.0*Math.Atan(num*srat(en.e*Math.Sin(slp.phi), -0.5*en.e))-Proj.HALFPI;
				if(Math.Abs(elp.phi-slp.phi)<TOL14) break;
				slp.phi=elp.phi;
			}

			// convergence failed
			if(i==0) Proj.pj_ctx_set_errno(ctx, -17);
			return elp;
		}
	}
}
