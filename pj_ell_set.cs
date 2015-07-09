using System;
using System.Collections.Generic;

// set ellipsoid parameters a and es

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// initialize geographic shape parameters
		internal static bool pj_ell_set(projCtx ctx, List<string> pl, out double a, out double es)
		{
			const double SIXTH=0.1666666666666666667; // 1/6
			const double RA4=0.04722222222222222222; // 17/360
			const double RA6=0.02215608465608465608; // 67/3024
			const double RV4=0.06944444444444444444; // 5/72
			const double RV6=0.04243827160493827160; // 55/1296

			// clear any previous error
			pj_ctx_set_errno(ctx, 0);

			// check for varying forms of ellipsoid input
			a=es=0.0;
			// R takes precedence
			if(pj_param_t(ctx, pl, "R")) a=pj_param_d(ctx, pl, "R");
			else
			{ // probable elliptical figure
				bool removeLast2=false;
				try
				{
					// check if ellps present and temporarily append its values to pl
					string name=pj_param_s(ctx, pl, "ellps");
					if(!string.IsNullOrEmpty(name))
					{
						string major, ell;
						if(GetEllipsoid(name, out major, out ell)==null) { pj_ctx_set_errno(ctx, -9); return true; }
						pl.Add(pj_mkparam(major));
						pl.Add(pj_mkparam(ell));
						removeLast2=true;
					}

					double b=0.0, e;
					a=pj_param_d(ctx, pl, "a");
					if(pj_param_t(ctx, pl, "es")) es=pj_param_d(ctx, pl, "es"); // eccentricity squared
					else if(pj_param_t(ctx, pl, "e"))
					{ // eccentricity
						e=pj_param_d(ctx, pl, "e");
						es=e*e;
					}
					else if(pj_param_t(ctx, pl, "rf"))
					{ // recip flattening
						es=pj_param_d(ctx, pl, "rf");
						if(es==0) { pj_ctx_set_errno(ctx, -10); return true; }
						es=1.0/es;
						es=es*(2.0-es);
					}
					else if(pj_param_t(ctx, pl, "f"))
					{ // flattening
						es=pj_param_d(ctx, pl, "f");
						es=es*(2.0-es);
					}
					else if(pj_param_t(ctx, pl, "b"))
					{ // minor axis
						b=pj_param_d(ctx, pl, "b");
						es=1.0-(b*b)/(a*a);
					}
					// else es == 0.0 and sphere of radius a

					if(b==0) b=a*Math.Sqrt(1.0-es);

					// following options turn ellipsoid into equivalent sphere
					if(pj_param_b(ctx, pl, "R_A"))
					{ // sphere--area of ellipsoid
						a*=1.0-es*(SIXTH+es*(RA4+es*RA6));
						es=0.0;
					}
					else if(pj_param_b(ctx, pl, "R_V"))
					{ // sphere--vol. of ellipsoid
						a*=1.0-es*(SIXTH+es*(RV4+es*RV6));
						es=0.0;
					}
					else if(pj_param_b(ctx, pl, "R_a"))
					{ // sphere--arithmetic mean
						a=0.5*(a+b);
						es=0.0;
					}
					else if(pj_param_b(ctx, pl, "R_g"))
					{ // sphere--geometric mean
						a=Math.Sqrt(a*b);
						es=0.0;
					}
					else if(pj_param_b(ctx, pl, "R_h"))
					{ // sphere--harmonic mean
						a=2.0*a*b/(a+b);
						es=0.0;
					}
					else
					{
						bool i=pj_param_t(ctx, pl, "R_lat_a"); // sphere--arith.
						if(i||pj_param_t(ctx, pl, "R_lat_g"))
						{ // or geom. mean at latitude
							double tmp=Math.Sin(pj_param_r(ctx, pl, i?"R_lat_a":"R_lat_g"));
							if(Math.Abs(tmp)>HALFPI) { pj_ctx_set_errno(ctx, -11); return true; }
							tmp=1.0-es*tmp*tmp;
							a*=i?0.5*(1.0-es+tmp)/(tmp*Math.Sqrt(tmp)):Math.Sqrt(1.0-es)/tmp;
							es=0.0;
						}
					}
				}
				finally
				{
					if(removeLast2) pl.RemoveRange(pl.Count-2, 2); // clean up temporary extension of list
				}

				if(ctx.last_errno!=0) return true;
			}

			// some remaining checks
			if(es<0.0) { pj_ctx_set_errno(ctx, -12); return true; }
			if(a<=0.0) { pj_ctx_set_errno(ctx, -13); return true; }

			return false;
		}
	}
}
