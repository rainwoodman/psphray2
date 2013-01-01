import ConfigParser
import argparse
import numpy

parser = argparse.ArgumentParser()
parser.add_argument("paramfile")
parser.add_argument("-q", '--cubic', action='store_true', default=False)

def parseargs(parser):
  args = parser.parse_args()
  config = ConfigParser.ConfigParser()
  config.read(args.paramfile)

  args.datadir = config.get("IO", "datadir")
  args.BoxSize = config.getfloat("IC", "BoxSize")
  args.a = config.getfloat("IC", "a")
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
  for name, value in config.items("Levels"):
     n, s = [float(t) for t in value.replace(';', ',').split(',')]
     Levels.append(n)
     args.Scale[n] = s
     args.Meta[n] = read_meta(n, args.datadir)

  args.Levels = numpy.array(Levels, dtype='i8')
  args.Levels.sort()

  return args

def read_meta(Nmesh, datadir):
    meta = {'Offset':{}, 'Size':{}}
    metalines = file('%s/filters-%d' % (datadir, Nmesh)).read()
    exec(metalines, meta)
    return meta


args = parseargs(parser)

def setup_level(Nmesh):
  if args.Scale[Nmesh] == 0: 
    return setup_baselevel(Nmesh)
  else:
    return setup_innerlevel(Nmesh)
  
def readblock(Nmesh, block, dtype):
  meta = args.Meta[Nmesh]
  NTask = meta['NTask']
  ZeroPadding = int(numpy.ceil(numpy.log10(NTask+0.1)))
  content = []
  for fid in range(NTask):
    content.append(numpy.fromfile('%s/%s-%d.%0*d' %(
                        args.datadir, 
                        block, Nmesh, ZeroPadding, fid), dtype))
  return numpy.concatenate(content)

def setup_innerlevel(Nmesh):
  index = readblock(Nmesh, 'index', 'i8')
  r = readblock(Nmesh, 'region', 'u1')
  dispx = readblock(Nmesh, 'dispx', 'f4')
  dispy = readblock(Nmesh, 'dispy', 'f4')
  dispz = readblock(Nmesh, 'dispz', 'f4')
  meta = args.Meta[Nmesh]
  sep = args.BoxSize / Nmesh
  scale = args.Scale[Nmesh]

  alldata = []
  for i, region in enumerate(args.Regions):
    mask = r == i
    data = numpy.empty(mask.sum(), dtype=[
           ('gps', ('i4', 3)), 
           ('disp', ('f4', 3))])
    data['gps'][:, 0], data['gps'][:, 1], data['gps'][:, 2] = numpy.unravel_index(index[mask], meta['Size'][i])
    data['gps'] += meta['Offset'][i]
    data['disp'][:, 0] = dispx[mask]
    data['disp'][:, 1] = dispy[mask]
    data['disp'][:, 2] = dispz[mask]
    alldata.append(data)
  data = numpy.concatenate(alldata)
  # now gonna remove the overlaps.
  numpy.remainder(data['gps'], Nmesh, data['gps'])
  index = numpy.ravel_multi_index(data['gps'].T, (Nmesh, Nmesh, Nmesh))
  u, ui = numpy.unique(index, return_index=True)
  return data[ui]

def setup_baselevel(Nmesh):
  data = numpy.empty((Nmesh, Nmesh, Nmesh), dtype=[
           ('gps', ('i4', 3)), 
           ('disp', ('f4', 3))])

  grid1d = numpy.arange(Nmesh);

  data['gps'][..., 2] = grid1d[None, None, :]
  data['gps'][..., 1] = grid1d[None, :, None]
  data['gps'][..., 0] = grid1d[:, None, None]

  data = data.reshape(-1)
  for i, block in enumerate(['dispx', 'dispy', 'dispz']):
    content = readblock(Nmesh, block, 'f4')
    data['disp'][..., i] = content

#  data['pos'] += data['vel']
#  data['vel'] *= meta['vfact']
#  numpy.remainder(data['pos'], args.BoxSize, data['pos'])
  return data

def dig(data, Nmesh, scale, Align=None):
  """ dig a hole on Nmesh points data, for NmeshDig
      returns the mask of the hole,
      align to Align, if omitted use Nmesh
  """
  if Align is None: Align = Nmesh
  sep = args.BoxSize / Align
  assert Nmesh % Align == 0
  mask2 = numpy.zeros(len(data), dtype='?')
  for region in args.Regions:
    gridcenter = region['center'] / sep
    dis = data['gps'] // (Nmesh // Align) - gridcenter
    # use the if any part of the mesh overlaps the sphere,
    dis[dis < - 0.5] += 1
    dis *= sep
    dis += args.BoxSize * 0.5
    numpy.remainder(dis, args.BoxSize, dis)
    dis -= args.BoxSize * 0.5
    dis /= (0.5 * region['size'] * scale)
    if not args.cubic:
      dis **= 2
      dis = dis.sum(axis=-1)
      mask2 |= (dis < 1.0)
    else:
      numpy.abs(dis, dis)
      mask2 |= (dis < 1.0).all(axis=-1)
  return mask2

def data2par(data, Nmesh):
  result = numpy.empty(shape=len(data), dtype=[
            ('pos', ('f4', 3)), ('vel', ('f4', 3)), ('Nmesh', 'i4')])
  result['pos'] = data['gps'] * args.BoxSize / Nmesh
  result['vel'] = args.Meta[Nmesh]['vfact'] * data['disp']
  result['Nmesh'] = Nmesh
  return result
i = 0
NmeshNext = args.Levels[0]
dataNext = setup_level(NmeshNext)
alldata = []
while i < len(args.Levels) - 1:
  Nmesh, data = NmeshNext, dataNext
  NmeshNext = args.Levels[i + 1]

  holemask = dig(data, Nmesh, scale=args.Scale[NmeshNext])
  print Nmesh, holemask.sum()
  alldata.append(data2par(data[~holemask], Nmesh))
  assert NmeshNext % Nmesh == 0

  dataNext = setup_level(NmeshNext)
  fillmask = dig(dataNext, NmeshNext, scale=args.Scale[NmeshNext], Align=Nmesh)
  dataNext = dataNext[fillmask]
  i = i + 1

Nmesh, data = NmeshNext, dataNext
alldata.append(data2par(data, Nmesh))


