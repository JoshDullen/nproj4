using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Gridshift
{
	public static partial class Grid
	{
		// Determine nad table correction value
		public static LP nad_intr(LP t, CTABLE ct)
		{
			LP val, frct;
			ILP indx;
			double m00, m10, m01, m11;
			int f00, f10, f01, f11;
			int index;
			int @in;

			t.lam/=ct.del.lam;
			t.phi/=ct.del.phi;
			indx.lam=(int)Math.Floor(t.lam);
			indx.phi=(int)Math.Floor(t.phi);
			frct.lam=t.lam-indx.lam;
			frct.phi=t.phi-indx.phi;
			val.lam=val.phi=Libc.HUGE_VAL;

			if(indx.lam<0)
			{
				if(indx.lam==-1&&frct.lam>0.99999999999)
				{
					indx.lam++;
					frct.lam=0.0;
				}
				else return val;
			}
			else
			{
				@in=indx.lam+1;
				if(@in>=ct.lim.lam)
				{
					if(@in==ct.lim.lam&&frct.lam<1e-11)
					{
						indx.lam--;
						frct.lam=1.0;
					}
					else return val;
				}
			}

			if(indx.phi<0)
			{
				if(indx.phi==-1&&frct.phi>0.99999999999)
				{
					indx.phi++;
					frct.phi=0.0;
				}
				else return val;
			}
			else
			{
				@in=indx.phi+1;
				if(@in>=ct.lim.phi)
				{
					if(@in==ct.lim.phi&&frct.phi<1e-11)
					{
						indx.phi--;
						frct.phi=1.0;
					}
					else return val;
				}
			}

			index=indx.phi*ct.lim.lam+indx.lam;
			f00=index++;
			f10=index;
			index+=ct.lim.lam;
			f11=index--;
			f01=index;
			m11=m10=frct.lam;
			m00=m01=1.0-frct.lam;
			m11*=frct.phi;
			m01*=frct.phi;
			frct.phi=1.0-frct.phi;
			m00*=frct.phi;
			m10*=frct.phi;

			val.lam=m00*ct.cvs[f00].lam+m10*ct.cvs[f10].lam+m01*ct.cvs[f01].lam+m11*ct.cvs[f11].lam;
			val.phi=m00*ct.cvs[f00].phi+m10*ct.cvs[f10].phi+m01*ct.cvs[f01].phi+m11*ct.cvs[f11].phi;

			return val;
		}
	}
}
