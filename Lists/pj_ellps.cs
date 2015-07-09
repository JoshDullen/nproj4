namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// definition of standard geoids
		public static string GetEllipsoid(string name, out string major, out string ell)
		{
			switch(name)
			{
				case "MERIT": major="a=6378137.0"; ell="rf=298.257"; return "MERIT 1983";
				case "SGS85": major="a=6378136.0"; ell="rf=298.257"; return "Soviet Geodetic System 85";
				case "GRS80": major="a=6378137.0"; ell="rf=298.257222101"; return "GRS 1980(IUGG, 1980)";
				case "IAU76": major="a=6378140.0"; ell="rf=298.257"; return "IAU 1976";
				case "airy": major="a=6377563.396"; ell="b=6356256.910"; return "Airy 1830";
				case "APL4.9": major="a=6378137.0"; ell="rf=298.25"; return "Appl. Physics. 1965";
				case "NWL9D": major="a=6378145.0"; ell="rf=298.25"; return "Naval Weapons Lab., 1965";
				case "mod_airy": major="a=6377340.189"; ell="b=6356034.446"; return "Modified Airy";
				case "andrae": major="a=6377104.43"; ell="rf=300.0"; return "Andrae 1876 (Den., Iclnd.)";
				case "aust_SA": major="a=6378160.0"; ell="rf=298.25"; return "Australian Natl & S. Amer. 1969";
				case "GRS67": major="a=6378160.0"; ell="rf=298.2471674270"; return "GRS 67(IUGG 1967)";
				case "bessel": major="a=6377397.155"; ell="rf=299.1528128"; return "Bessel 1841";
				case "bess_nam": major="a=6377483.865"; ell="rf=299.1528128"; return "Bessel 1841 (Namibia)";
				case "clrk66": major="a=6378206.4"; ell="b=6356583.8"; return "Clarke 1866";
				case "clrk80": major="a=6378249.145"; ell="rf=293.4663"; return "Clarke 1880 mod.";
				case "clrk80ign": major="a=6378249.2"; ell="rf=293.4660212936269"; return "Clarke 1880 (IGN)";
				case "CPM": major="a=6375738.7"; ell="rf=334.29"; return "Comm. des Poids et Mesures 1799";
				case "delmbr": major="a=6376428.0"; ell="rf=311.5"; return "Delambre 1810 (Belgium)";
				case "engelis": major="a=6378136.05"; ell="rf=298.2566"; return "Engelis 1985";
				case "evrst30": major="a=6377276.345"; ell="rf=300.8017"; return "Everest 1830";
				case "evrst48": major="a=6377304.063"; ell="rf=300.8017"; return "Everest 1948";
				case "evrst56": major="a=6377301.243"; ell="rf=300.8017"; return "Everest 1956";
				case "evrst69": major="a=6377295.664"; ell="rf=300.8017"; return "Everest 1969";
				case "evrstSS": major="a=6377298.556"; ell="rf=300.8017"; return "Everest (Sabah & Sarawak)";
				case "fschr60": major="a=6378166.0"; ell="rf=298.3"; return "Fischer (Mercury Datum) 1960";
				case "fschr60m": major="a=6378155.0"; ell="rf=298.3"; return "Modified Fischer 1960";
				case "fschr68": major="a=6378150.0"; ell="rf=298.3"; return "Fischer 1968";
				case "helmert": major="a=6378200.0"; ell="rf=298.3"; return "Helmert 1906";
				case "hough": major="a=6378270.0"; ell="rf=297.0"; return "Hough";
				case "intl": major="a=6378388.0"; ell="rf=297.0"; return "International 1909 (Hayford)";
				case "krass": major="a=6378245.0"; ell="rf=298.3"; return "Krassovsky, 1942";
				case "kaula": major="a=6378163.0"; ell="rf=298.24"; return "Kaula 1961";
				case "lerch": major="a=6378139.0"; ell="rf=298.257"; return "Lerch 1979";
				case "mprts": major="a=6397300.0"; ell="rf=191.0"; return "Maupertius 1738";
				case "new_intl": major="a=6378157.5"; ell="b=6356772.2"; return "New International 1967";
				case "plessis": major="a=6376523.0"; ell="b=6355863.0"; return "Plessis 1817 (France)";
				case "SEasia": major="a=6378155.0"; ell="b=6356773.3205"; return "Southeast Asia";
				case "walbeck": major="a=6376896.0"; ell="b=6355834.8467"; return "Walbeck";
				case "WGS60": major="a=6378165.0"; ell="rf=298.3"; return "WGS 60";
				case "WGS66": major="a=6378145.0"; ell="rf=298.25"; return "WGS 66";
				case "WGS72": major="a=6378135.0"; ell="rf=298.26"; return "WGS 72";
				case "WGS84": major="a=6378137.0"; ell="rf=298.257223563"; return "WGS 84";
				case "sphere": major="a=6370997.0"; ell="b=6370997.0"; return "Normal Sphere (r=6370997)";
			}

			major=ell=null;
			return null;
		}
	}
}
