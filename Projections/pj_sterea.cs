// libproj -- library of cartographic projections
//
// Copyright (c) 2003 Gerald I. Evenden
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
using System;
using Free.Ports.Proj4.Gauss;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_sterea : PJ
	{
		protected double phic0, cosc0, sinc0, R2;
		protected GAUSS en;

		public override string Name { get { return "sterea"; } }
		public override string DescriptionName { get { return "Oblique Stereographic Alternative"; } }
		public override string DescriptionType { get { return "Azi, Sph&Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			lp=GAUSS.pj_gauss(ctx, lp, en);
			double sinc=Math.Sin(lp.phi);
			double cosc=Math.Cos(lp.phi);
			double cosl=Math.Cos(lp.lam);
			double k=k0*R2/(1.0+sinc0*sinc+cosc0*cosc*cosl);
			xy.x=k*cosc*Math.Sin(lp.lam);
			xy.y=k*(cosc0*sinc-sinc0*cosc*cosl);

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.x/=k0;
			xy.y/=k0;
			double rho=Libc.hypot(xy.x, xy.y);
			if(rho!=0)
			{
				double c=2.0*Math.Atan2(rho, R2);
				double sinc=Math.Sin(c);
				double cosc=Math.Cos(c);
				lp.phi=Math.Asin(cosc*sinc0+xy.y*sinc*cosc0/rho);
				lp.lam=Math.Atan2(xy.x*sinc, rho*cosc0*cosc-xy.y*sinc0*sinc);
			}
			else
			{
				lp.phi=phic0;
				lp.lam=0.0;
			}

			return GAUSS.pj_inv_gauss(ctx, lp, en);
		}

		public override PJ Init()
		{
			double R;

			en=GAUSS.pj_gauss_ini(e, phi0, out phic0, out R);
			if(en==null) return null;

			sinc0=Math.Sin(phic0);
			cosc0=Math.Cos(phic0);
			R2=2.0*R;
			inv=e_inverse;
			fwd=e_forward;

			return this;
		}
	}
}
