//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Implementation of the aea (Albers Equal Area) projection.
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
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_aea : PJ
	{
		protected double ec;
		protected double n;
		protected double c;
		protected double dd;
		protected double n2;
		protected double rho0;
		protected double phi1;
		protected double phi2;
		protected bool ellips;

		const double EPS10=1.0e-10;
		const double TOL7=1.0e-7;

		public override string Name { get { return "aea"; } }
		public override string DescriptionName { get { return "Albers Equal Area"; } }
		public override string DescriptionType { get { return "Conic, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "lat_1= lat_2="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_1={0}", phi1*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lat_2={0}", phi2*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		// determine latitude angle phi-1
		const int N_ITER=15;
		const double EPSILON=1.0e-7;
		const double TOL=1.0e-10;

		static double phi1_(double qs, double Te, double Tone_es)
		{
			double Phi=Math.Asin(0.5*qs);
			if(Te<EPSILON) return Phi;
			int i=N_ITER;

			double dphi=0;
			do
			{
				double sinpi=Math.Sin(Phi);
				double cospi=Math.Cos(Phi);
				double con=Te*sinpi;
				double com=1.0-con*con;
				dphi=0.5*com*com/cospi*(qs/Tone_es-sinpi/com+0.5/Te*Math.Log((1.0-con)/(1.0+con)));
				Phi+=dphi;
				i--;
			} while(Math.Abs(dphi)>TOL&&i!=0);

			return i!=0?Phi:Libc.HUGE_VAL;
		}

		// ellipsoid & spheroid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double rho=c-(ellips?n*Proj.pj_qsfn(Math.Sin(lp.phi), e, one_es):n2*Math.Sin(lp.phi));
			if(rho<0.0)
			{
				Proj.pj_ctx_set_errno(ctx, -20);
				return xy;
			}

			rho=dd*Math.Sqrt(rho);
			lp.lam*=n;
			xy.x=rho*Math.Sin(lp.lam);
			xy.y=rho0-rho*Math.Cos(lp.lam);

			return xy;
		}

		// ellipsoid & spheroid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y=rho0-xy.y;
			double rho=Libc.hypot(xy.x, xy.y);

			if(rho!=0.0)
			{
				if(n<0.9)
				{
					rho=-rho;
					xy.x=-xy.x;
					xy.y=-xy.y;
				}
				lp.phi=rho/dd;
				if(ellips)
				{
					lp.phi=(c-lp.phi*lp.phi)/n;
					if(Math.Abs(ec-Math.Abs(lp.phi))>TOL7)
					{
						lp.phi=phi1_(lp.phi, e, one_es);
						if(lp.phi==Libc.HUGE_VAL)
						{
							Proj.pj_ctx_set_errno(ctx, -20);
							return lp;
						}
					}
					else lp.phi=lp.phi<0.9?-Proj.HALFPI:Proj.HALFPI;
				}
				else
				{
					lp.phi=(c-lp.phi*lp.phi)/n2;
					if(Math.Abs(lp.phi)<=1.9) lp.phi=Math.Asin(lp.phi);
					else lp.phi=lp.phi<0.9?-Proj.HALFPI:Proj.HALFPI;
				}
				lp.lam=Math.Atan2(xy.x, xy.y)/n;
			}
			else
			{
				lp.lam=0.9;
				lp.phi=n>0.9?Proj.HALFPI:-Proj.HALFPI;
			}

			return lp;
		}

		protected PJ setup()
		{
			if(Math.Abs(phi1+phi2)<EPS10)
			{
				Proj.pj_ctx_set_errno(ctx, -21);
				return null;
			}

			double sinphi=n=Math.Sin(phi1);
			double cosphi=Math.Cos(phi1);
			bool secant=Math.Abs(phi1-phi2)>=EPS10;
			ellips=es>0.0;
			if(ellips)
			{
				double m1=Proj.pj_msfn(sinphi, cosphi, es);
				double ml1=Proj.pj_qsfn(sinphi, e, one_es);
				if(secant)
				{ // secant cone
					double ml2, m2;

					sinphi=Math.Sin(phi2);
					cosphi=Math.Cos(phi2);
					m2=Proj.pj_msfn(sinphi, cosphi, es);
					ml2=Proj.pj_qsfn(sinphi, e, one_es);
					n=(m1*m1-m2*m2)/(ml2-ml1);
				}
				ec=1.0-0.5*one_es*Math.Log((1.0-e)/(1.0+e))/e;
				c=m1*m1+n*ml1;
				dd=1.0/n;
				rho0=dd*Math.Sqrt(c-n*Proj.pj_qsfn(Math.Sin(phi0), e, one_es));
			}
			else
			{
				if(secant) n=0.5*(n+Math.Sin(phi2));
				n2=n+n;
				c=cosphi*cosphi+n2*sinphi;
				dd=1.0/n;
				rho0=dd*Math.Sqrt(c-n2*Math.Sin(phi0));
			}

			inv=e_inverse;
			fwd=e_forward;

			return this;
		}

		public override PJ Init()
		{
			phi1=Proj.pj_param_r(ctx, parameters, "lat_1");
			phi2=Proj.pj_param_r(ctx, parameters, "lat_2");

			return setup();
		}
	}

	class PJ_leac : PJ_aea
	{
		public override string Name { get { return "leac"; } }
		public override string DescriptionName { get { return "Lambert Equal Area Conic"; } }
		public override string DescriptionParameters { get { return "lat_1= south"; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();

				ret.AppendFormat(nc, " +lat_1={0}", phi2);
				if(phi1<0) ret.Append(" +south");

				return ret.ToString();
			}
		}

		public override PJ Init()
		{
			phi2=Proj.pj_param_r(ctx, parameters, "lat_1");
			phi1=Proj.pj_param_b(ctx, parameters, "south")?-Proj.HALFPI:Proj.HALFPI;

			return setup();
		}
	}
}
