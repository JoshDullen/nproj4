using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_chamb : PJ
	{
		protected struct VECT { public double r, Az; }

		protected struct C
		{ // control point data
			public double phi, lam;
			public double cosphi, sinphi;
			public VECT v;
			public XY p;
		}

		protected C[] c=new C[3];
		protected XY p;
		protected double beta_0, beta_1, beta_2;

		public override string Name { get { return "chamb"; } }
		public override string DescriptionName { get { return "Chamberlin Trimetric"; } }
		public override string DescriptionType { get { return "Misc Sph, no inv."; } }
		public override string DescriptionParameters { get { return "lat_1= lon_1= lat_2= lon_2= lat_3= lon_3="; } }
		public override bool Invertible { get { return false; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_1={0}", c[0].phi*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lon_1={0}", Proj.adjlon(c[0].lam+lam0)*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lat_2={0}", c[1].phi*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lon_2={0}", Proj.adjlon(c[1].lam+lam0)*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lat_3={0}", c[2].phi*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lon_3={0}", Proj.adjlon(c[2].lam+lam0)*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		const double THIRD=0.333333333333333333;
		const double TOL9=1.0e-9;

		// distance and azimuth from point 1 to point 2
		static VECT vect(projCtx ctx, double dphi, double c1, double s1, double c2, double s2, double dlam)
		{
			VECT v;

			double cdl=Math.Cos(dlam);
			if(Math.Abs(dphi)>1.0||Math.Abs(dlam)>1.0) v.r=Proj.aacos(ctx, s1*s2+c1*c2*cdl);
			else
			{ // more accurate for smaller distances
				double dp=Math.Sin(0.5*dphi);
				double dl=Math.Sin(0.5*dlam);
				v.r=2.0*Proj.aasin(ctx, Math.Sqrt(dp*dp+c1*c2*dl*dl));
			}
			if(Math.Abs(v.r)>TOL9) v.Az=Math.Atan2(c2*Math.Sin(dlam), c1*s2-s1*c2*cdl);
			else v.r=v.Az=0.0;

			return v;
		}

		// law of cosines
		static double lc(projCtx ctx, double b, double c, double a)
		{
			return Proj.aacos(ctx, 0.5*(b*b+c*c-a*a)/(b*c));
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			VECT[] v=new VECT[3];

			double sinphi=Math.Sin(lp.phi);
			double cosphi=Math.Cos(lp.phi);

			int i=0;
			for(; i<3; ++i)
			{ // dist/azimiths from control
				v[i]=vect(ctx, lp.phi-c[i].phi, c[i].cosphi, c[i].sinphi, cosphi, sinphi, lp.lam-c[i].lam);
				if(v[i].r==0) break;
				v[i].Az=Proj.adjlon(v[i].Az-c[i].v.Az);
			}

			if(i<3) xy=c[i].p; // current point at control point
			else
			{ // point mean of intersepts
				xy=p;
				for(i=0; i<3; ++i)
				{
					int j=(i==2)?0:i+1;
					double a=lc(ctx, c[i].v.r, v[i].r, v[j].r);
					if(v[i].Az<0.0) a=-a;

					if(i==0)
					{ // coord comp unique to each arc
						xy.x+=v[i].r*Math.Cos(a);
						xy.y-=v[i].r*Math.Sin(a);
					}
					else if(i==1)
					{
						a=beta_1-a;
						xy.x-=v[i].r*Math.Cos(a);
						xy.y-=v[i].r*Math.Sin(a);
					}
					else
					{
						a=beta_2-a;
						xy.x+=v[i].r*Math.Cos(a);
						xy.y+=v[i].r*Math.Sin(a);
					}
				}
				xy.x*=THIRD; // mean of arc intercepts
				xy.y*=THIRD;
			}

			return xy;
		}

		public override PJ Init()
		{
			for(int i=0; i<3; i++)
			{ // get control point locations
				string line="lat_"+(i+1).ToString();
				c[i].phi=Proj.pj_param_r(ctx, parameters, line);
				line="lon_"+(i+1).ToString();
				c[i].lam=Proj.pj_param_r(ctx, parameters, line);
				c[i].lam=Proj.adjlon(c[i].lam-lam0);
				c[i].cosphi=Math.Cos(c[i].phi);
				c[i].sinphi=Math.Sin(c[i].phi);
			}

			for(int i=0; i<3; i++)
			{ // inter ctl pt. distances and azimuths
				int j=i==2?0:i+1;
				c[i].v=vect(ctx, c[j].phi-c[i].phi, c[i].cosphi, c[i].sinphi, c[j].cosphi, c[j].sinphi, c[j].lam-c[i].lam);
				if(c[i].v.r==0) { Proj.pj_ctx_set_errno(ctx, -25); return null; }
				// co-linearity problem ignored for now
			}

			beta_0=lc(ctx, c[0].v.r, c[2].v.r, c[1].v.r);
			beta_1=lc(ctx, c[0].v.r, c[1].v.r, c[2].v.r);
			beta_2=Proj.PI-beta_0;
			p.y=2.0*(c[0].p.y=c[1].p.y=c[2].v.r*Math.Sin(beta_0));
			c[2].p.y=0.0;
			c[0].p.x=-(c[1].p.x=0.5*c[0].v.r);
			p.x=c[2].p.x=c[0].p.x+c[2].v.r*Math.Cos(beta_0);
			es=0.0; fwd=s_forward;

			return this;
		}
	}
}
