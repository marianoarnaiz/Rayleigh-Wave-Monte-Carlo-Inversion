function cm_data=tempe(m)
cm = [[0.027451, 0.117647, 0.27451],
[0.027451, 0.184314, 0.419608],
[0.0313725, 0.321569, 0.611765],
[0.129412, 0.443137, 0.709804],
[0.258824, 0.572549, 0.780392],
[0.352941, 0.627451, 0.803922],
[0.470588, 0.74902, 0.839216],
[0.666667, 0.862745, 0.901961],
[0.858824, 0.960784, 1],
[0.941176, 0.988235, 1],
[1, 0.941176, 0.960784],
[1, 0.878431, 0.878431],
[0.988235, 0.733333, 0.666667],
[0.988235, 0.572549, 0.447059],
[0.984314, 0.415686, 0.290196],
[0.941176, 0.235294, 0.168627],
[0.8, 0.0941176, 0.117647],
[0.65098, 0.0588235, 0.0784314],
[0.470588, 0.0392157, 0.0588235],
[0.372549, 0, 0]];

if nargin < 1
    cm_data = cm;
else
    hsv=rgb2hsv(cm);
    cm_data=interp1(linspace(0,1,size(cm,1)),hsv,linspace(0,1,m));
    cm_data=hsv2rgb(cm_data);
  
end
end