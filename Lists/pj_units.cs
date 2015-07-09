namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// definition of standard cartesian units

		// returns the multiplier to convert named units to meters
		// may be expressed by either a simple floating point constant or a
		// numerator/denomenator values (e.g. 1/1000)
		public static double GetUnitFactor(string name)
		{
			switch(name)
			{
				case "km": return 1000.0;
				case "m": return 1.0;
				case "dm": return 0.1;
				case "cm": return 0.01;
				case "mm": return 0.001;
				case "kmi": return 1852.0;
				case "in": return 0.0254;
				case "ft": return 0.3048;
				case "yd": return 0.9144;
				case "mi": return 1609.344;
				case "fath": return 1.8288;
				case "ch": return 20.1168;
				case "link": return 0.201168;
				case "us-in": return 1.0/39.37;
				case "us-ft": return 0.304800609601219;
				case "us-yd": return 0.914401828803658;
				case "us-ch": return 20.11684023368047;
				case "us-mi": return 1609.347218694437;
				case "ind-yd": return 0.91439523;
				case "ind-ft": return 0.30479841;
				case "ind-ch": return 20.11669506;
				default: return double.NaN;
			}
		}

		public static string GetUnitName(string name)
		{
			switch(name)
			{
				case "km": return "Kilometer";
				case "m": return "Meter";
				case "dm": return "Decimeter";
				case "cm": return "Centimeter";
				case "mm": return "Millimeter";
				case "kmi": return "International Nautical Mile";
				case "in": return "International Inch";
				case "ft": return "International Foot";
				case "yd": return "International Yard";
				case "mi": return "International Statute Mile";
				case "fath": return "International Fathom";
				case "ch": return "International Chain";
				case "link": return "International Link";
				case "us-in": return "U.S. Surveyor's Inch";
				case "us-ft": return "U.S. Surveyor's Foot";
				case "us-yd": return "U.S. Surveyor's Yard";
				case "us-ch": return "U.S. Surveyor's Chain";
				case "us-mi": return "U.S. Surveyor's Statute Mile";
				case "ind-yd": return "Indian Yard";
				case "ind-ft": return "Indian Foot";
				case "ind-ch": return "Indian Chain";
				default: return null;
			}
		}
	}
}
