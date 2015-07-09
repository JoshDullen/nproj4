//*****************************************************************************
//
// Project: PROJ.4
// Purpose: Implementation of the nzmg (New Zealand Map Grid) projection.
//			Very loosely based upon DMA code by Bradford W. Drew
// Author:	Gerald Evenden
//
//*****************************************************************************
// Copyright (c) 1995, Gerald Evenden
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
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//*****************************************************************************

using System;
using Free.Ports.Proj4.ComplexPoly;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_nzmg : PJ
	{
		public override string Name { get { return "nzmg"; } }
		public override string DescriptionName { get { return "New Zealand Map Grid"; } }
		public override string DescriptionType { get { return "fixed Earth"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double EPSLN=1.0e-10;
		const double SEC5_TO_RAD=0.4848136811095359935899141023;
		const double RAD_TO_SEC5=2.062648062470963551564733573;

		static readonly COMPLEX[] bf=new COMPLEX[]
		{
			new COMPLEX(0.7557853228, 0.0),
			new COMPLEX(0.249204646, 0.003371507),
			new COMPLEX(-0.001541739, 0.041058560),
			new COMPLEX(-0.10162907, 0.01727609),
			new COMPLEX(-0.26623489, -0.36249218),
			new COMPLEX(-0.6870983, -1.1651967)
		};

		static readonly double[] tphi=new double[] { 1.5627014243, 0.5185406398, -0.03333098, -0.1052906, -0.0368594, 0.007317, 0.01220, 0.00394, -0.0013 };
		static readonly double[] tpsi=new double[] { 0.6399175073, -0.1358797613, 0.063294409, -0.02526853, 0.0117879, -0.0055161, 0.0026906, -0.001333, 0.00067, -0.00034 };

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			COMPLEX p;

			lp.phi=(lp.phi-phi0)*RAD_TO_SEC5;

			int i=tpsi.Length-1;
			int C=i;
			for(p.r=tpsi[C--]; i>0; i--) p.r=tpsi[C--]+lp.phi*p.r;

			p.r*=lp.phi;
			p.i=lp.lam;
			p=p.pj_zpoly1(bf, bf.Length-1);
			xy.x=p.i;
			xy.y=p.r;

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			COMPLEX p;
			p.r=xy.y;
			p.i=xy.x;

			int nn;
			for(nn=20; nn>0; nn--)
			{
				COMPLEX fp;
				COMPLEX f=p.pj_zpolyd1(bf, bf.Length-1, out fp);
				f.r-=xy.y;
				f.i-=xy.x;
				double den=fp.r*fp.r+fp.i*fp.i;
				COMPLEX dp;
				dp.r=-(f.r*fp.r+f.i*fp.i)/den;
				dp.i=-(f.i*fp.r-f.r*fp.i)/den;
				p.r+=dp.r;
				p.i+=dp.i;
				if((Math.Abs(dp.r)+Math.Abs(dp.i))<=EPSLN) break;
			}

			if(nn!=0)
			{
				lp.lam=p.i;
				int i=tphi.Length-1;
				int C=i;
				for(lp.phi=tphi[C--]; i>0; i--) lp.phi=tphi[C--]+p.r*lp.phi;

				lp.phi=phi0+p.r*lp.phi*SEC5_TO_RAD;
			}
			else lp.lam=lp.phi=Libc.HUGE_VAL;

			return lp;
		}

		public override PJ Init()
		{
			// force to International major axis
			a=6378388.0;
			ra=1.0/a;
			lam0=Proj.DEG_TO_RAD*173.0;
			phi0=Proj.DEG_TO_RAD*-41.0;
			x0=2510000.0;
			y0=6023150.0;
			inv=e_inverse;
			fwd=e_forward;

			return this;
		}
	}
}
