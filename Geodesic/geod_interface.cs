using System;

namespace Free.Ports.Proj4.Geodesic
{
	public partial class Geod
	{
		// was GEODESIC.A ... GEODESIC.DIST
		public double geod_a=0, geod_f=0, lam1=0, phi1=0, al12=0, lam2=0, phi2=0, al21=0, geod_S=0;

		public geod_geodesic GlobalGeodesic;
		public geod_geodesicline GlobalGeodesicLine;
		public int n_alpha, n_S;
		public double to_meter, fr_meter, del_alpha;

		public void Ini()
		{
			GlobalGeodesic=new geod_geodesic(geod_a, geod_f);
		}

		public void Pre()
		{
			double degree=Math.PI/180;
			double lat1=phi1/degree, lon1=lam1/degree, azi1=al12/degree;

			GlobalGeodesicLine=new geod_geodesicline(GlobalGeodesic, lat1, lon1, azi1, GEOD.NONE);
		}

		public void For()
		{
			double degree=Math.PI/180, s12=geod_S;

			double lat2, lon2, azi2;
			GlobalGeodesicLine.geod_position(s12, out lat2, out lon2, out azi2);

			azi2+=azi2>=0?-180:180; // Compute back azimuth
			phi2=lat2*degree;
			lam2=lon2*degree;
			al21=azi2*degree;
		}

		public void Inv()
		{
			double degree=Math.PI/180;
			double lat1=phi1/degree, lon1=lam1/degree, lat2=phi2/degree, lon2=lam2/degree;
			
			double azi1, azi2, s12;
			GlobalGeodesic.geod_inverse(lat1, lon1, lat2, lon2, out s12, out azi1, out azi2);

			azi2+=azi2>=0?-180:180; // Compute back azimuth
			al12=azi1*degree;
			al21=azi2*degree;
			geod_S=s12;
		}
	}
}
