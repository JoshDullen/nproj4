//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Implementation of the airy (Airy) projection.
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

namespace Free.Ports.Proj4.Projections
{
	class PJ_airy : PJ
	{
		protected double p_halfpi;
		protected double sinph0;
		protected double cosph0;
		protected double Cb;
		protected airy_mode mode;
		protected bool no_cut; // do not cut at hemisphere limit

		public override string Name { get { return "airy"; } }
		public override string DescriptionName { get { return "Airy"; } }
		public override string DescriptionType { get { return "Misc Sph, no inv."; } }
		public override string DescriptionParameters { get { return "no_cut lat_b="; } }
		public override bool Invertible { get { return false; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();

				if(no_cut) ret.Append(" +no_cut");

				double x=Proj.HALFPI/2;
				if(Cb<=-0.5+EPS) x=0;
				else if(Math.Abs(Cb)<EPS) x=Proj.HALFPI;
				else
				{
					double diff=0;
					do
					{
						double sinx=Math.Sin(x);
						double cosx=Math.Cos(x);
						double logcosx=Math.Log(cosx);
						diff=(cosx*cosx*sinx*logcosx-Cb*sinx*sinx*sinx)/(2*cosx*logcosx+cosx*sinx*sinx);
						x=x+diff;
					} while(Math.Abs(diff)>EPS);
				}

				ret.AppendFormat(nc, " +lat_b={0}", (Proj.HALFPI-2*x)*Proj.RAD_TO_DEG);

				return ret.ToString();
			}
		}

		const double EPS=1.0e-10;

		protected enum airy_mode
		{
			N_POLE=0,
			S_POLE=1,
			EQUIT=2,
			OBLIQ=3
		}

		// spheroid
		XY s_forward_airy(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double sinlam=Math.Sin(lp.lam);
			double coslam=Math.Cos(lp.lam);

			switch(mode)
			{
				case airy_mode.EQUIT:
				case airy_mode.OBLIQ:
					{
						double sinphi=Math.Sin(lp.phi);
						double cosphi=Math.Cos(lp.phi);
						double cosz=cosphi*coslam;
						if(mode==airy_mode.OBLIQ) cosz=sinph0*sinphi+cosph0*cosz;
						if(!no_cut&&cosz<-EPS)
						{
							Proj.pj_ctx_set_errno(ctx, -20);
							return xy;
						}

						double s=1.0-cosz;
						double Krho;
						if(Math.Abs(s)>EPS)
						{
							double t=0.5*(1.0+cosz);
							Krho=-Math.Log(t)/s-Cb/t;
						}
						else Krho=0.5-Cb;

						xy.x=Krho*cosphi*sinlam;
						if(mode==airy_mode.OBLIQ) xy.y=Krho*(cosph0*sinphi-sinph0*cosphi*coslam);
						else xy.y=Krho*sinphi;
					}
					break;
				case airy_mode.S_POLE:
				case airy_mode.N_POLE:
					{
						lp.phi=Math.Abs(p_halfpi-lp.phi);
						if(!no_cut&&(lp.phi-EPS)>Proj.HALFPI)
						{
							Proj.pj_ctx_set_errno(ctx, -20);
							return xy;
						}

						lp.phi*=0.5;
						if(lp.phi>EPS)
						{
							double t=Math.Tan(lp.phi);
							double Krho=-2.0*(Math.Log(Math.Cos(lp.phi))/t+t*Cb);
							xy.x=Krho*sinlam;
							xy.y=Krho*coslam;
							if(mode==airy_mode.N_POLE) xy.y=-xy.y;
						}
						else xy.x=xy.y=0.0;
					}
					break;
			}
			return xy;
		}

		public override PJ Init()
		{
			no_cut=Proj.pj_param_b(ctx, parameters, "no_cut");
			double beta=0.5*(Proj.HALFPI-Proj.pj_param_r(ctx, parameters, "lat_b"));
			if(Math.Abs(beta)<EPS) Cb=-0.5;
			else
			{
				Cb=1.0/Math.Tan(beta);
				Cb*=Cb*Math.Log(Math.Cos(beta));
			}

			if(Math.Abs(Math.Abs(phi0)-Proj.HALFPI)<EPS)
			{
				if(phi0<0.0)
				{
					p_halfpi=-Proj.HALFPI;
					mode=airy_mode.S_POLE;
				}
				else
				{
					p_halfpi=Proj.HALFPI;
					mode=airy_mode.N_POLE;
				}
			}
			else
			{
				if(Math.Abs(phi0)<EPS) mode=airy_mode.EQUIT;
				else
				{
					mode=airy_mode.OBLIQ;
					sinph0=Math.Sin(phi0);
					cosph0=Math.Cos(phi0);
				}
			}
			fwd=s_forward_airy;
			es=0.0;

			return this;
		}
	}
}
