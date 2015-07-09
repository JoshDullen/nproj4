//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Implementation of the aitoff (Aitoff) and wintri (Winkel Tripel)
//			projections.
// Author:	Gerald Evenden
//
//*****************************************************************************
// Copyright (c) 1995, Gerald Evenden
// Copyright (c) 2008-2015 by the Authors
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
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_aitoff : PJ
	{
		protected double cosphi1;
		protected bool mode;

		public override string Name { get { return "aitoff"; } }
		public override string DescriptionName { get { return "Aitoff"; } }
		public override string DescriptionType { get { return "Misc Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double c=0.5*lp.lam;
			double d=Math.Acos(Math.Cos(lp.phi)*Math.Cos(c));

			if(d!=0) // basic Aitoff
			{
				xy.y=1.0/Math.Sin(d);
				xy.x=2.0*d*Math.Cos(lp.phi)*Math.Sin(c)*xy.y;
				xy.y*=d*Math.Sin(lp.phi);
			}
			else xy.x=xy.y=0;

			if(mode) // Winkel Tripel
			{
				xy.x=(xy.x+lp.lam*cosphi1)*0.5;
				xy.y=(xy.y+lp.phi)*0.5;
			}

			return xy;
		}

		//***********************************************************************************
		//
		// Inverse functions added by Drazen Tutic and Lovro Gradiser based on paper:
		//
		// I.Özbug Biklirici and Cengizhan Ipbüker. A General Algorithm for the Inverse
		// Transformation of Map Projections Using Jacobian Matrices. In Proceedings of the
		// Third International Symposium Mathematical & Computational Applications,
		// pages 175-182, Turkey, September 2002.
		//
		// Expected accuracy is defined by EPSILON = 1e-12. Should be appropriate for
		// most applications of Aitoff and Winkel Tripel projections.
		//
		// Longitudes of 180W and 180E can be mixed in solution obtained.
		//
		// Inverse for Aitoff projection in poles is undefined, longitude value of 0 is assumed.
		//
		// Contact : dtutic@geof.hr
		// Date: 2015-02-16
		//
		//***********************************************************************************
		const double EPSILON=1e-12;
		const double HALFPI=Math.PI/2;
		const int MAXITER=10, MAXROUND=20;

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			if(Math.Abs(xy.x)<EPSILON&&Math.Abs(xy.y)<EPSILON) return lp;

			int iter, round=0;
			double x, y, dp, dl;

			// intial values for Newton-Raphson method
			lp.lam=xy.x;
			lp.phi=xy.y;

			do
			{
				iter=0;
				do
				{
					double sl=Math.Sin(lp.lam*0.5);
					double cl=Math.Cos(lp.lam*0.5);
					double sp=Math.Sin(lp.phi);
					double cp=Math.Cos(lp.phi);

					double d=cp*cl;
					double c=1.0-d*d;
					d=Math.Acos(d)/Math.Pow(c, 1.5);

					double f1=2.0*d*c*cp*sl;
					double f2=d*c*sp;
					double f1p=2.0*(sl*cl*sp*cp/c-d*sp*sl);
					double f1l=cp*cp*sl*sl/c+d*cp*cl*sp*sp;
					double f2p=sp*sp*cl/c+d*sl*sl*cp;
					double f2l=0.5*(sp*cp*sl/c-d*sp*cp*cp*sl*cl);

					if(mode) // Winkel Tripel
					{
						f1=0.5*(f1+lp.lam*cosphi1);
						f2=0.5*(f2+lp.phi);
						f1p*=0.5;
						f1l=0.5*(f1l+cosphi1);
						f2p=0.5*(f2p+1.0);
						f2l*=0.5;
					}

					f1-=xy.x;
					f2-=xy.y;

					dp=f1p*f2l-f2p*f1l;
					dl=(f2*f1p-f1*f2p)/dp;
					dp=(f1*f2l-f2*f1l)/dp;

					// set to interval [-π, π]
					while(dl>Math.PI) dl-=Math.PI;
					while(dl<-Math.PI) dl+=Math.PI;

					lp.lam-=dl;
					lp.phi-=dp;
				} while((Math.Abs(dp)>EPSILON||Math.Abs(dl)>EPSILON)&&iter++<MAXITER);

				// correct if symmetrical solution for Aitoff
				if(lp.phi>HALFPI) lp.phi-=2.0*(lp.phi-HALFPI);
				if(lp.phi<-HALFPI) lp.phi-=2.0*(lp.phi+HALFPI);

				// if pole in Aitoff, return longitude of 0
				if((Math.Abs(Math.Abs(lp.phi)-HALFPI)<EPSILON)&&!mode) lp.lam=0.0;

				// calculate x, y coordinates with solution obtained
				double C=0.5*lp.lam;
				double D=Math.Acos(Math.Cos(lp.phi)*Math.Cos(C));

				if(D!=0) // Aitoff
				{
					y=1.0/Math.Sin(D);
					x=2.0*D*Math.Cos(lp.phi)*Math.Sin(C)*y;
					y*=D*Math.Sin(lp.phi);
				}
				else x=y=0;

				if(mode) // Winkel Tripel
				{
					x=(x+lp.lam*cosphi1)*0.5;
					y=(y+lp.phi)*0.5;
				}

				// if too far from given values of x, y, repeat with better approximation of phi, lam
			} while((Math.Abs(xy.x-x)>EPSILON||Math.Abs(xy.y-y)>EPSILON)&&round++<MAXROUND);

			if(iter==MAXITER&&round==MAXROUND) Proj.pj_log(ctx, PJ_LOG.ERROR, "Warning: Accuracy of 1e-12 not reached. Last increments: dlat={0} and dlon={0}", dp, dl);

			return lp;
		}

		protected PJ setup()
		{
			inv=s_inverse;
			fwd=s_forward;
			es=0.0;
			return this;
		}

		public override PJ Init()
		{
			mode=false;

			return setup();
		}
	}

	class PJ_wintri : PJ_aitoff
	{
		public override string Name { get { return "wintri"; } }
		public override string DescriptionName { get { return "Winkel Tripel"; } }
		public override string DescriptionType { get { return "Misc Sph"; } }
		public override string DescriptionParameters { get { return "lat_1="; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_1={0}", Math.Acos(cosphi1)*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		public override PJ Init()
		{
			mode=true;
			if(Proj.pj_param_t(ctx, parameters, "lat_1"))
			{
				cosphi1=Math.Cos(Proj.pj_param_r(ctx, parameters, "lat_1"));
				if(cosphi1==0.0)
				{
					Proj.pj_ctx_set_errno(ctx, -22);
					return null;
				}
			}
			else cosphi1=0.636619772367581343; // 50d28' or Math.Acos(2/pi)

			return setup();
		}
	}
}
