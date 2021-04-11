import math


class BoltGroup:
    def __init__(self):
        #self.coords = [(-150,150,1),(-50,150,2), (50,150,3), (150,150,4),(-150,50,5),(-50,50,6), (50,50,7), (150,50,8),(-150,-50,9),(-50,-50,10), (50,-50,11), (150,-50,12),(-150,-150,13),(-50,-150,14), (50,-150,15), (150,-150,16)]
        #9 x M20
        #self.coords = [(-100,100,1),(0,100,2),(100,100,3),(-100,0,1),(0,0,2),(100,0,3),(-100,-100,1),(0,-100,2),(100,-100,3)]
        #25 x M7
        #self.coords = [(-70,70,2),(-35,70,3),(0,70,4),(35,70,5),(70,70,6),(-70,35,2),(-35,35,3),(0,35,4),(35,35,5),(70,35,6),(-70,0,2),(-35,0,3),(0,0,4),(35,0,5),(70,0,6),(-70,-35,2),(-35,-35,3),(0,-35,4),(35,-35,5),(70,-35,6),(-70,-70,2),(-35,-70,3),(0,-70,4),(35,-70,5),(70,-70,6)]
        #20 x M7
        self.coords = [(0, 30, 1), (0,-30,2)]
        #54 x m7
        #self.coords = [(-85.5,140,1),(-52.5,140,2),(-17.5,140,3),(17.5,140,4),(52.5,140,5),(85.5,140,6),(-85.5,105,1),(-52.5,105,2),(-17.5,105,3),(17.5,105,4),(52.5,105,5),(85.5,105,6),(-85.5,70,1),(-52.5,70,2),(-17.5,70,3),(17.5,70,4),(52.5,70,5),(85.5,70,6),(-85.5,35,1),(-52.5,35,2),(-17.5,35,3),(17.5,35,4),(52.5,35,5),(85.5,35,6),(-85.5,0,1),(-52.5,0,2),(-17.5,0,3),(17.5,0,4),(52.5,0,5),(85.5,0,6),(-85.5,-35,1),(-52.5,-35,2),(-17.5,-35,3),(17.5,-35,4),(52.5,-35,5),(85.5,-35,6),(-85.5,-70,1),(-52.5,-70,2),(-17.5,-70,3),(17.5,-70,4),(52.5,-70,5),(85.5,-70,6),(-85.5,-105,1),(-52.5,-105,2),(-17.5,-105,3),(17.5,-105,4),(52.5,-105,5),(85.5,-105,6),(-85.5,-140,1),(-52.5,-140,2),(-17.5,-140,3),(17.5,-140,4),(52.5,-140,5),(85.5,-140,6)]
        self.grain_direction = 90
        self.bolt_num = len(self.coords)
        self.bolt_size = 12
        self.forcey = -6.5
        self.forcex = 0
        self.moment = 1.2
        self.k1 = 0.77
        self.phi = 0.8
        self.plate = 1
        self.width = 130
        self.k16 = 1
        self.k17 = 1
        self.timber_density = 500
        self.plate_thickness = 12
        self.jointgroup = "JD4"
        self.parawidth = 90
        self.perpwidth = 70
        self.CLT = False
        if self.plate == 1:
            self.effectivewidth = (self.width - self.plate_thickness)/2
            self.effectiveparawidth = (self.parawidth - self.plate_thickness/2)/2
            self.effectiveperpwidth = self.perpwidth/2
        else:
            self.effectivewidth = self.width - 2*self.plate_thickness
            self.effectiveperpwidth = self.perpwidth - 2*self.plate_thickness
            self.effectiveparawidth = self.parawidth

        if self.grain_direction == 90:
            self.forcey, self.forcex =  self.forcex, self.forcey


    def centroid(self):
        SumX = 0
        SumY = 0
        for i in range(self.bolt_num):
            SumX += self.coords[i][0]
            SumY += self.coords[i][1]

        X = SumX / self.bolt_num
        Y = SumY / self.bolt_num
        return (X,Y)

    def calculateforces(self):
        centroid = self.centroid()
        cx = centroid[0]
        cy = centroid[1]
        Icx = 0
        Icy = 0
        forcelist = []

        for i in range(self.bolt_num):
            Rcx = self.coords[i][0] - cx
            Rcy = self.coords[i][1] - cy
            Icx += Rcx**2
            Icy += Rcy**2

        polarmoment = (Icx + Icy)
        
        for i in range(self.bolt_num): 
            Rcx = self.coords[i][0] - cx
            Rcy = self.coords[i][1] - cy
            ForceXBolt = round(self.forcex / self.bolt_num - ((self.moment * Rcy * 1000) / polarmoment),2);
            ForceYBolt = round(self.forcey / self.bolt_num + ((self.moment * Rcx * 1000) / polarmoment),2);
            forcelist.append((ForceXBolt,ForceYBolt,self.coords[i][2]))

        return forcelist

    def findQkl(self,width):
        fcjdict = {
            "J1": 55.5,
            "J2": 44,
            "J3": 35.5,
            "J4": 28,
            "J5": 22,
            "J6": 18,
            "JD1": 69,
            "JD2": 55.5,
            "JD3": 44,
            "JD4": 35.5,
            "JD5": 28,
            "JD6": 22
        }

        fcj = fcjdict[self.jointgroup]

        Qkl1 = width * fcj * self.bolt_size / 2;

        if "1" in self.jointgroup:
            Qkl2 = 1.65 * fcj * self.bolt_size**2

        elif "2" in self.jointgroup:
            Qkl2 = 1.75 * fcj * self.bolt_size**2

        elif "3" in self.jointgroup or "4" in self.jointgroup: 
            Qkl2 = 2 * fcj * self.bolt_size**2

        elif "5" in self.jointgroup: 
            Qkl2 = 2.2 * fcj * self.bolt_size**2

        else:
            Qkl2 = 2.4 * fcj * self.bolt_size**2

        Qkl = min(Qkl1, Qkl2)

        return Qkl/1000

    def findQkp(self,width):
        fpjdict = {

            "J1": 22,
            "J2": 17.5, 
            "J3": 11, 
            "J4": 7.1,
            "J5": 4.7,
            "J6": 2.4,
            "JD1": 29.5, 
            "JD2": 22.5,
            "JD3": 17,
            "JD4": 12.5, 
            "JD5": 9,
            "JD6": 6.1 
        }

        fpj = fpjdict[self.jointgroup]
        Qkp1 = self.width * fpj * self.bolt_size / 2

        if "1" in self.jointgroup:
            Qkp2 = 10 * fpj * math.sqrt(self.bolt_size**3)
        elif "2" in self.jointgroup:
            Qkp2 = 12 * fpj * math.sqrt(self.bolt_size**3)
        elif "3" in self.jointgroup:
            Qkp2 = 15 * fpj * math.sqrt(self.bolt_size**3)
        elif "4" in self.jointgroup:
            Qkp2 = 17 * fpj * math.sqrt(self.bolt_size**3)
        elif "5" in self.jointgroup:
            Qkp2 = 19 * fpj * math.sqrt(self.bolt_size**3)
        
        else:
            Qkp2 = 22 * fpj * math.sqrt(self.bolt_size**3)

        Qkp = min(Qkp1, Qkp2)
        return Qkp/1000

    def CheckCapacity(self):
        forcelist = self.calculateforces()
        capacitylist = []
        if self.CLT == False:
            Qkp = self.findQkp(self.effectivewidth)
            Qkl = self.findQkl(self.effectivewidth)
        else:
            Qkp_temp = self.findQkp(self.effectivewidth)
            Qkl_temp = self.findQkl(self.effectivewidth)
            pararatio = self.effectiveparawidth/self.effectivewidth
            perpratio = self.effectiveperpwidth/self.effectivewidth
            Qkl = Qkp_temp*perpratio + Qkl_temp*pararatio
            Qkp = Qkp_temp*pararatio + Qkl_temp*perpratio

        for boltforce in forcelist:
            if boltforce[0] == 0:
                angle = 90
            elif boltforce[1] == 0:
                angle = 0
            elif self.grain_direction == 0:
                angle = math.degrees(math.atan(abs(boltforce[1]/boltforce[0])))
            else:
                angle = math.degrees(math.atan(abs(boltforce[0]/boltforce[1])))
            force = math.sqrt(boltforce[0]**2 + boltforce[1]**2)
            if boltforce[1] == 0:
                Qk = Qkl
            elif boltforce[0] == 0:
                Qk = Qkp
            else:
                Qk = (Qkl * Qkp) / ((Qkl * math.sin(math.radians(angle))**2) + (Qkp * math.cos(math.radians(angle))**2))


            boltcapacity = round((force / (2 * self.phi * self.k1 * self.k16 * self.k17 * Qk)*100),2)

            capacitylist.append((boltcapacity, boltforce[2]))

        return capacitylist

    def findDowelStiffness(self):
        # kNm/rad
        k_ser = (2*self.timber_density)**(1.5)*self.bolt_size/23
        k_u = (2/3)*k_ser
        summation = 0
        centroid = self.centroid()
        cx = centroid[0]
        cy = centroid[1]
        for i in range(len(self.coords)):
            Rcx = self.coords[i][0] - cx
            Rcy = self.coords[i][1] - cy
            summation += math.sqrt(Rcx**2+Rcy**2)**2
        stiffness_rotation = k_ser*summation/10**6

        return stiffness_rotation

bg = BoltGroup()
print(bg.CheckCapacity())
print("stiffness_rotation is " + str(bg.findDowelStiffness()))