import ConfigParser
import argparse
import numpy
import os.path
import StringIO

# this requires numpy > 1.6.0

def parseargs():
  parser = argparse.ArgumentParser(description="Assembling IC files from dumps")
  parser.add_argument("paramfile", 
               help="the same paramfile fed to ngenic.grid.")
  parser.add_argument("format", choices=['gadget', 'ramses'], 
                 help="Format of the output. \
                       gadget : to write a gadget IC, \
                       ramses : to write a ramses IC, only the first zoom region is written due to ramses limitation")
  parser.add_argument("-N", "--num-particles", 
                 help="Number of DM particles per file. \
                       For levels with Gas the total number in a file will be twice of this",
                      default=1024 * 1024 * 4, dest='Npar', type=int);
  parser.add_argument("--iddtype",  choices=['uint32', 'uint64'],
                 help="bit width of the id field",
                 default=numpy.dtype(numpy.uint64), dest='iddtype', type=numpy.dtype)
  
  parser.add_argument("--idscheme",  choices=['seq', 'mesh'],
                 help="Format of the id. 'serial' is to sequentially \
                       assign particle ids. \
                       for 'mesh', the 64 bit id = 4 bits empty,\
                                            1 bit  gas(1)/dm(0), \
                                            7 bit level index, \
                               52 bits serial mesh index for that level ",
                       default='mesh', dest='idscheme')
                       
  parser.add_argument("-o", "--prefix", 
                 help="Prefix of the IC files. \
                       Do not forget to put a filename base. \
                       gadget : IC.%%d, one file per zoom level. \
                       ramses : IC-%%d-%%s, several files per zoom level",
                                    default="IC")
  parser.add_argument("-q", '--cubic', 
                 help="Use an axis-aligned box for the zoom region, instead of spheres", 
                     action='store_true', default=False)
  parser.add_argument("-B", "--base", nargs='?',
                 help="treat the base level like the finest level, \
                       and skip all other levels. \
                       In other words, write the entire base level mesh, \
                       into particle type specified in the finest level, \
                       generating gas if the finest level says so. \
                       accepts --base=N where the base will be downsampled to Nmesh/N. \
                       ", 
                 action="store", default=False, const=1)

  parser.add_argument('--no-displacement', dest='nodisp', 
                  help="Do not displace the position of particles. \
                        Physics will be corrupted. Do not use", 
                        action='store_true', default=False)
  args = parser.parse_args()
  config = ConfigParser.ConfigParser()
  config.readfp(StringIO.StringIO(
"""[Cosmology]
h = 0.72
H = 0.1
G = 43007.1
C = 3e5
OmegaB = 0.044
OmegaM = 0.26
OmegaL = 0.74
[Units]
UnitLength_in_cm = 3.085678e21
UnitMass_in_g = 1.989e43
UnitVelocity_in_cm_per_s = 1e5
"""))
  str = file(args.paramfile).read().replace(';', ',').replace('#', ';')
  config.readfp(StringIO.StringIO(str))
  args.datadir = config.get("IO", "datadir")
  args.BoxSize = config.getfloat("IC", "BoxSize")
  args.a = config.getfloat("IC", "a")
  args.OmegaM = config.getfloat("Cosmology", "OmegaM")
  args.OmegaB = config.getfloat("Cosmology", "OmegaB")
  args.OmegaL = config.getfloat("Cosmology", "OmegaL")
  args.h = config.getfloat("Cosmology", "h")
  args.UnitLength_in_cm = config.getfloat("Unit", "UnitLength_in_cm")
  args.UnitVelocity_in_cm_per_s = config.getfloat("Unit", "UnitVelocity_in_cm_per_s")
  args.UnitMass_in_g = config.getfloat("Unit", "UnitMass_in_g")
  args.UnitTime_in_s = args.UnitLength_in_cm / args.UnitVelocity_in_cm_per_s;
  args.G = config.getfloat("Cosmology", "G")
  args.H = config.getfloat("Cosmology", "H")

  Regions = []
  for name, value in config.items("Regions"):
    all = [float(t) for t in value.replace(';', ',').split(',')]
    if len(all) == 4:
      x, y, z, s = all
      sx = s
      sy = s
      sz = s
    elif len(all) == 6:
      x, y, z, sx, sy, sz = all
  
    Regions.append(((x, y, z), (sx, sy, sz)))
  args.Regions = numpy.array(Regions, dtype=[('center', ('f8', 3)), ('size', ('f8', 3))])
  
  Levels = []
  
  args.Scale = {}
  args.Meta = {}
  args.DMptype = {}
  args.MakeGas = {}
  args.DownSample = {}
  for name, value in config.items("Levels"):
    sp = value.replace(';', ',').split(',')
    good = False
    if len(sp) < 2 or len(sp) > 4:
      raise 'paramfile format error, need 2/3/4 entries per Level [Levels]'

    if len(sp) >= 2:
      n, s = int(sp[0]), float(sp[1])
      # the default is 
      p = -1
      g = None
      d = 1
    if len(sp) >= 3:
      if 'g' in sp[2]:
        p = int(sp[2].replace('g', ''))
        g = True
      else:
        p = int(sp[2])
        g = False
    if len(sp) == 4:
      d = int(sp[3])

    Levels.append(n)
    args.Scale[n] = s
    args.DMptype[n] = p
    args.MakeGas[n] = g
    args.DownSample[n] = d
    args.Meta[n] = read_meta(n, args.datadir)

  args.Levels = numpy.array(Levels, dtype='i8')
  args.Levels.sort()

  assert args.Scale[args.Levels[0]] == 0.0

  for Nmesh in args.Levels:
    if args.DMptype[Nmesh] == -1:
      if Nmesh == args.Levels.max():
        args.DMptype[Nmesh] = 1
        args.MakeGas[Nmesh] = True
      else:
        args.DMptype[Nmesh] = 2
        args.MakeGas[Nmesh] = False

  if args.base: 
    # set up the base Nmesh to coarses Nmesh / args.base,
    # and rely on readblock to scale to the new Nmesh / args.base
    # file downsample
    first = args.Levels[0]
    # must divide the base level
    assert first % args.base == 0
    Base = first / args.base
    last = args.Levels[-1]
    args.DMptype[Base] = args.DMptype[last]
    args.MakeGas[Base] = args.MakeGas[last]
    args.Levels = numpy.array([Base], dtype='i8')
    args.Meta[Base] = args.Meta[first]
  return args

def read_meta(Nmesh, datadir):
    meta = {'Offset':{}, 'Size':{}}
    metalines = file('%s/meta-%d' % (datadir, Nmesh)).read()
    exec(metalines, meta)
    return meta

def scale(Nmesh, src, ratio):
  Nmeshsrc = Nmesh * ratio
  # must be covering the full space
  assert Nmeshsrc ** 3 == len(src)
  src = src.reshape(Nmesh, ratio, Nmesh, ratio, Nmeshsrc, ratio)
  src = src.sum(axis=(5, 3, 1)).ravel()
  src /= (ratio ** 3)
  return src
def readblock(Nmesh, block, dtype, args):
  if args.base:
    Nmesh = Nmesh * args.base
    # doesn't make sense no region and index block
    # on the base level.
    # and scaling them is nonsense.
    assert block != 'region' 
    assert block != 'index'
      

  header = {}
  exec(
    file('%s/%s-%d.header' % (args.datadir, block, Nmesh)).read(), 
    header)
  NTask = header['NTask']
  DownSample = header['DownSample']
  # otherwise getting a wrong file.
  assert DownSample == args.DownSample[Nmesh]

  width1 = len(str(DownSample))
  width2 = len(str(NTask))
  content = []

  f = 0
  for i in range(DownSample):
    for fid in range(NTask):
      if DownSample > 1:
        fname = '%s/%s-%d.%0*d-%0*d' %(
                 args.datadir, 
                 block, Nmesh, width1, i, width2, fid)
      else:
        fname = '%s/%s-%d.%0*d' %(
                 args.datadir, 
                   block, Nmesh, width2, fid)
      if os.path.exists(fname):
        content.append(numpy.fromfile(fname, dtype))
        f = f + 1
      #skip the tasks that do not dump any files
  print Nmesh, block, f, 'files read'
  content = numpy.concatenate(content)
  if len(content) == 0:
    raise Exception("the Nmesh = %d has no content(too small), has to terminate"
                % Nmesh)
  if args.base:
    content = scale(Nmesh, content, args.base)

  return content

def ramses(args):
  def write_ramses(header, Nmesh, name, buffer):
    with file('%s%d-%s' % (args.prefix, Nmesh, name), 'w') as icfile:
      buffer = buffer.reshape(header['np'])
      writerecord(icfile, header)
      for row in buffer.transpose((2, 0, 1)):
        writerecord(icfile, row.ravel('F'))

  RHEADER = numpy.dtype([
   ('np', ('i4', 3)), ('dx', 'f4'), 
   ('xo', ('f4', 3)), ('astart', 'f4'), 
   ('omegam', 'f4'), ('omegav', 'f4'), ('h0', 'f4')])
  fac = 1.0 / args.h * args.UnitLength_in_cm / 3.08567758e24
  boxsize = args.BoxSize * fac
  for i, Nmesh in enumerate(args.Levels):
    header = numpy.zeros(None, RHEADER)
    if args.Scale[Nmesh] == 0:
      header['np'][:] = Nmesh
      header['xo'][:] = 0
      mask = None
    else:
      header['np'][:] = args.Meta[Nmesh]['Size'][0]
      header['xo'][:] = header['dx'] * args.Meta[Nmesh]['Offset'][0]
      r = readblock(Nmesh, 'region', 'u4', args)
      mask = r == 0
    header['dx'] = boxsize / Nmesh
    header['astart'] = args.a
    header['omegam'] = args.OmegaM
    header['omegav'] = args.OmegaL
    header['h0'] = args.h * 100.0
    print Nmesh, header

    delta = readblock(Nmesh, 'delta', 'f4', args)
    if mask is not None: delta = delta[mask]
    write_ramses(header, Nmesh, 'deltab', delta)
    delta = None
    for i, block in enumerate(['dispx', 'dispy', 'dispz']):
      disp = readblock(Nmesh, block, 'f4', args)
      if mask is not None: disp = disp[mask]
      disp *= fac
      write_ramses(header, Nmesh, block, disp)
      disp /= args.Meta[Nmesh]['vfact']
      write_ramses(header, Nmesh, ['velcx', 'velcy', 'velcz'][i], disp)

def main(args):
  if args.format == 'gadget':
    from gadgetmodule import gadget
    gadget(args)
  elif args.format == 'ramses':
    ramses(args)

if __name__ == '__main__':
  main(parseargs())
