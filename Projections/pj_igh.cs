using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_igh : PJ
	{
		protected PJ[] pj=new PJ[12];
		protected double dy0;

		public override string Name { get { return "igh"; } }
		public override string DescriptionName { get { return "Interrupted Goode Homolosine"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double d4044118=(40+44/60.0+11.8/3600.0)*Proj.DEG_TO_RAD; // 40d 44' 11.8" [degrees]

		const double d10=10*Proj.DEG_TO_RAD;
		const double d20=20*Proj.DEG_TO_RAD;
		const double d30=30*Proj.DEG_TO_RAD;
		const double d40=40*Proj.DEG_TO_RAD;
		const double d50=50*Proj.DEG_TO_RAD;
		const double d60=60*Proj.DEG_TO_RAD;
		const double d80=80*Proj.DEG_TO_RAD;
		const double d90=90*Proj.DEG_TO_RAD;
		const double d100=100*Proj.DEG_TO_RAD;
		const double d140=140*Proj.DEG_TO_RAD;
		const double d160=160*Proj.DEG_TO_RAD;
		const double d180=180*Proj.DEG_TO_RAD;

		const double EPSLN=1.0e-10; // allow a little 'slack' on zone edge positions

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			int z=0;
			if(lp.phi>=d4044118) z=(lp.lam<=-d40?1:2); // 1|2
			else if(lp.phi>=0) z=(lp.lam<=-d40?3:4); // 3|4
			else if(lp.phi>=-d4044118)
			{ // 5|6|7|8
				if(lp.lam<=-d100) z=5; // 5
				else if(lp.lam<=-d20) z=6; // 6
				else if(lp.lam<=d80) z=7; // 7
				else z=8; // 8
			}
			else
			{ // 9|10|11|12
				if(lp.lam<=-d100) z=9; // 9
				else if(lp.lam<=-d20) z=10; // 10
				else if(lp.lam<=d80) z=11; // 11
				else z=12; // 12
			}

			lp.lam-=pj[z-1].lam0;
			xy=pj[z-1].fwd(lp);
			xy.x+=pj[z-1].x0;
			xy.y+=pj[z-1].y0;

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double y90=dy0+Math.Sqrt(2); // lt=90 corresponds to y=y0+sqrt(2)

			int z=0;
			if(xy.y>y90+EPSLN||xy.y<-y90+EPSLN) z=0;// 0
			else if(xy.y>=d4044118) z=(xy.x<=-d40?1:2); // 1|2
			else if(xy.y>=0) z=(xy.x<=-d40?3:4); // 3|4
			else if(xy.y>=-d4044118)
			{ // 5|6|7|8
				if(xy.x<=-d100) z=5; // 5
				else if(xy.x<=-d20) z=6; // 6
				else if(xy.x<=d80) z=7; // 7
				else z=8; // 8
			}
			else
			{ // 9|10|11|12
				if(xy.x<=-d100) z=9; // 9
				else if(xy.x<=-d20) z=10; // 10
				else if(xy.x<=d80) z=11; // 11
				else z=12; // 12
			}

			if(z!=0)
			{
				bool ok=false;

				xy.x-=pj[z-1].x0;
				xy.y-=pj[z-1].y0;
				lp=pj[z-1].inv(xy);
				lp.lam+=pj[z-1].lam0;

				switch(z)
				{
					case 1: ok=(lp.lam>=-d180-EPSLN&&lp.lam<=-d40+EPSLN)||((lp.lam>=-d40-EPSLN&&lp.lam<=-d10+EPSLN)&&
							(lp.phi>=d60-EPSLN&&lp.phi<=d90+EPSLN)); break;
					case 2: ok=(lp.lam>=-d40-EPSLN&&lp.lam<=d180+EPSLN)||((lp.lam>=-d180-EPSLN&&lp.lam<=-d160+EPSLN)&&
							(lp.phi>=d50-EPSLN&&lp.phi<=d90+EPSLN))||((lp.lam>=-d50-EPSLN&&lp.lam<=-d40+EPSLN)&&
							(lp.phi>=d60-EPSLN&&lp.phi<=d90+EPSLN)); break;
					case 3: ok=(lp.lam>=-d180-EPSLN&&lp.lam<=-d40+EPSLN); break;
					case 4: ok=(lp.lam>=-d40-EPSLN&&lp.lam<=d180+EPSLN); break;
					case 5: ok=(lp.lam>=-d180-EPSLN&&lp.lam<=-d100+EPSLN); break;
					case 6: ok=(lp.lam>=-d100-EPSLN&&lp.lam<=-d20+EPSLN); break;
					case 7: ok=(lp.lam>=-d20-EPSLN&&lp.lam<=d80+EPSLN); break;
					case 8: ok=(lp.lam>=d80-EPSLN&&lp.lam<=d180+EPSLN); break;
					case 9: ok=(lp.lam>=-d180-EPSLN&&lp.lam<=-d100+EPSLN); break;
					case 10: ok=(lp.lam>=-d100-EPSLN&&lp.lam<=-d20+EPSLN); break;
					case 11: ok=(lp.lam>=-d20-EPSLN&&lp.lam<=d80+EPSLN); break;
					case 12: ok=(lp.lam>=d80-EPSLN&&lp.lam<=d180+EPSLN); break;
				}

				z=(!ok?0:z); // projectable?
			}
			//if(z==0) { Proj.pj_ctx_set_errno(ctx, -15); return null; } // invalid x or y
			if(z==0) lp.lam=Libc.HUGE_VAL;
			if(z==0) lp.phi=Libc.HUGE_VAL;
			return lp;
		}

		//  Zones:
		//
		//    -180            -40                       180
		//      +--------------+-------------------------+    Zones 1,2,9,10,11 & 12:
		//      |1             |2                        |      Mollweide projection
		//      |              |                         |
		//      +--------------+-------------------------+    Zones 3,4,5,6,7 & 8:
		//      |3             |4                        |      Sinusoidal projection
		//      |              |                         |
		//    0 +-------+------+-+-----------+-----------+
		//      |5      |6       |7          |8          |
		//      |       |        |           |           |
		//      +-------+--------+-----------+-----------+
		//      |9      |10      |11         |12         |
		//      |       |        |           |           |
		//      +-------+--------+-----------+-----------+
		//    -180    -100      -20         80          180

		bool SETUP(int n, PJ proj, double x_0, double y_0, double lon_0)
		{
			if(proj==null) return false;
			pj[n-1]=proj.Init();
			if(pj[n-1]==null) return false;

			pj[n-1].x0=x_0;
			pj[n-1].y0=y_0;
			pj[n-1].lam0=lon_0;

			return true;
		}

		public override PJ Init()
		{
			LP lp=new LP();
			lp.lam=0;
			lp.phi=d4044118;

			XY xy1;
			XY xy3;

			try
			{
				// sinusoidal zones
				pj[3-1]=new PJ_sinu();
				pj[3-1].Init();

				if(!SETUP(3, new PJ_sinu(), -d100, 0, -d100)) return null;
				if(!SETUP(4, new PJ_sinu(), d30, 0, d30)) return null;
				if(!SETUP(5, new PJ_sinu(), -d160, 0, -d160)) return null;
				if(!SETUP(6, new PJ_sinu(), -d60, 0, -d60)) return null;
				if(!SETUP(7, new PJ_sinu(), d20, 0, d20)) return null;
				if(!SETUP(8, new PJ_sinu(), d140, 0, d140)) return null;

				// mollweide zones
				if(!SETUP(1, new PJ_moll(), -d100, 0, -d100)) return null;

				// y0 ?
				xy1=pj[0].fwd(lp); // zone 1
				xy3=pj[2].fwd(lp); // zone 3
				// y0 + xy1.y = xy3.y for lt = 40d44'11.8"
				dy0=xy3.y-xy1.y;

				pj[0].y0=dy0;

				// mollweide zones (cont'd)
				if(!SETUP(2, new PJ_moll(), d30, dy0, d30)) return null;
				if(!SETUP(9, new PJ_moll(), -d160, -dy0, -d160)) return null;
				if(!SETUP(10, new PJ_moll(), -d60, -dy0, -d60)) return null;
				if(!SETUP(11, new PJ_moll(), d20, -dy0, d20)) return null;
				if(!SETUP(12, new PJ_moll(), d140, -dy0, d140)) return null;
			}
			catch
			{
				return null;
			}

			inv=s_inverse;
			fwd=s_forward;
			es=0.0;

			return this;
		}
	}
}
