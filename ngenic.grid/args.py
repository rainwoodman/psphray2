import numpy
import os.path
import StringIO
import ConfigParser
import argparse

def parseargs():
  parser = argparse.ArgumentParser(description="Assembling IC files from dumps")
  parser.add_argument("paramfile", 
               help="the same paramfile fed to ngenic.grid.")
  parser.add_argument("action", choices=['powerspectrum', 'interpolate', 'gadget', 'ramses'], 
                 help="Format of the output. \
                       powerspectrum: write a pycamb powerspectrum file. \
                       interpolate: to fill the DownSample regions, \
                       gadget : to write a gadget IC, \
                       ramses : to write a ramses IC")
  parser.add_argument("-N", "--num-particles", 
                 help="Number of DM particles per file. \
                       For levels with Gas the total number \
                       in a file will be twice of this",
                      default=1024 * 1024 * 4, dest='Npar', type=int);

  parser.add_argument("--filter",  
                 help="order of the filter in interpolation",
                 default=4, dest='filter', type=int)
  parser.add_argument("--iddtype",  choices=['uint32', 'uint64'],
                 help="bit width of the id field",
                 default=numpy.dtype(numpy.uint64), dest='iddtype', type=numpy.dtype)
  
  parser.add_argument("-o", "--prefix", 
                 help="Prefix of the IC files. \
                       Do not forget to put a filename base. \
                       gadget : IC.%%d, one file per zoom level. \
                       ramses : IC-%%d-%%s, several files per zoom level",
                                    default="IC")
  parser.add_argument("-q", '--cubic', 
                 help="Use an axis-aligned box for the zoom region, instead of spheres", 
                     action='store_true', default=False)

  parser.add_argument('--no-displacement', dest='nodisp', 
                  help="Do not displace the position of particles. \
                        Physics will be corrupted. Do not use", 
                        action='store_true', default=False)
  args = parser.parse_args()
  config1 = ConfigParser.ConfigParser()
  ssss = file(args.paramfile).read().replace(';', ',').replace('#', ';')
  config1.readfp(StringIO.StringIO(ssss))
  datadir = config1.get("IO", "datadir")
  config = ConfigParser.ConfigParser()
  ssss = file(datadir + '/paramfile-used').read().replace(';', ',').replace('#', ';')
  # disable reading from -used
#  config.readfp(StringIO.StringIO(ssss))
  config.optionxform = str
  # paramfile-used contains defaults, yet
  # paramfile overrides paramfile-used,
  # sometimes we can pick datadir from an old run and just modify
  # the radius
  ssss = file(args.paramfile).read().replace(';', ',').replace('#', ';')
  config.readfp(StringIO.StringIO(ssss))
  args.datadir = config.get("IO", "datadir")
  args.BoxSize = config.getfloat("IC", "BoxSize")
  args.a = config.getfloat("IC", "a")
  args.NmeshPrimary = config.getfloat("IC", "NmeshPrimary")
  print args.NmeshPrimary
  args.OmegaM = config.getfloat("Cosmology", "OmegaM")
  args.OmegaB = config.getfloat("Cosmology", "OmegaB")
  args.OmegaL = config.getfloat("Cosmology", "OmegaL")
  args.Sigma8 = config.getfloat("Cosmology", "Sigma8")
  args.PrimordialIndex = config.getfloat("Cosmology", "PrimordialIndex")
  args.h = config.getfloat("Cosmology", "h")
  args.UnitLength_in_cm = config.getfloat("Unit", "UnitLength_in_cm")
  args.UnitVelocity_in_cm_per_s = config.getfloat("Unit", "UnitVelocity_in_cm_per_s")
  args.UnitMass_in_g = config.getfloat("Unit", "UnitMass_in_g")
  args.UnitTime_in_s = args.UnitLength_in_cm / args.UnitVelocity_in_cm_per_s;
  args.G = config.getfloat("Cosmology", "G")
  args.H = config.getfloat("Cosmology", "H")

  args.WhichSpectrum = config.getint("Cosmology", "WhichSpectrum");
  if args.WhichSpectrum == 2:
      args.PowerSpectrumFile = config.get("Cosmology", "PowerSpectrumFile")
  Regions = []
  for name, value in config.items("Regions"):
    all = [float(t) for t in value.split(',')]
    if len(all) == 4:
      x, y, z, s = all
      sx, sy, sz = s, s, s
    elif len(all) == 6:
      x, y, z, sx, sy, sz = all
  
    Regions.append(((x, y, z), (sx, sy, sz)))
  args.R = numpy.array(Regions, dtype=[('center', ('f8', 3)), ('size', ('f8', 3))])
  print args.R
  assert len(args.R) == 1 
  Levels = []
  
  for name, value in config.items("Levels"):
    sp = value.replace(';', ',').split(',')

    if len(sp) >= 4:
      # the default is 
      n, s, p, g, d = int(sp[0]), float(sp[1]), -1, False, 1
      if 'g' in sp[2]:
        p, g = int(sp[2].replace('g', '')), True
      else:
        p, g = int(sp[2]), False
      d = int(sp[3])
    else:
      raise Exception('need four entries per Level')
    Levels.append((n, s, p, g, d))
  
  args.L = numpy.array(Levels, dtype=[
           ('Nmesh', 'i4'),
           ('scale', 'f8'),
           ('ptype', 'i4'),
           ('makegas', '?'),
           ('downsample', 'i4')])
  args.L.sort(order='Nmesh')

  print args.L
  if numpy.sum(args.L['ptype'] == 1) > 1:
    print args.L
    raise Exception('More than 1 dm species assigned to ptype 1, gadget FOF will not'
                 + ' work properly')
  args.META = {}
  for i, l in enumerate(args.L):
    args.META[l['Nmesh']] = build_meta_block(args, i)

  return args

def build_meta_block(A, ilevel):
  l = A.L[ilevel]
  meta = {'Offset':{}, 'Size':{}}
  metalines = file('%s/delta-%d.header' % (A.datadir, l['Nmesh'])).read()
  exec(metalines, globals(), meta)
  meta['ilevel'] = ilevel
  meta['Offset'] = numpy.array(meta['Offset'][0], dtype='intp')
  meta['Size'] = numpy.array(meta['Size'][0], dtype='intp')
  meta['Nmesh'] = l['Nmesh']
  meta['NmeshPrimary'] = A.NmeshPrimary
  meta['Scale'] = l['scale']
  meta['datadir'] = A.datadir
  meta['ZoomCenter'] = A.R['center'][0]
  meta['ZoomRadius'] = A.R['size'][0] * 0.5 * l['scale']
  meta['ptype'] = l['ptype']
  meta['makegas'] = l['makegas']
  meta['A'] = A
  return meta
