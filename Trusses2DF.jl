#-------------------------------------------------------------------------------------------------------------------
# Consultas a cesar.lifonzo@unsch.edu.pe.
#--------------------------------------------------------------------------------------------------------------------

using ExcelReaders
using ExcelFiles, DataFrames
#using ExcelFiles, DataTables, IndexedTables, TimeSeries, Temporal
using CSV
#using RDatasets

#using Plots

#using PyPlot

# IMPORTING DATA
f = openxl("DATA.xlsx")
#data = readxlsheet("DATOS.xlsx", "Hoja1")
#data1 = readxl(f, "Nudos!B2:C7")
#data3 = DataFrame(load("DATOS.xlsx", "Nudos"))
#data4 = load("DATOS.xlsx", "Nudos")
nodes= readxlsheet(f, "Nodes")
nodes=nodes[2:end,:]
members= readxlsheet(f, "Members") #Start/End node of member
members=members[2:end,:]
properties= readxlsheet(f, "Properties")
properties=properties[2:end,:]
forces= readxlsheet(f, "Forces")
forces=forces[2:end,:]#Node Load Concenter
supported= readxlsheet(f, "Supported")
supported=supported[2:end,:]

nn=size(nodes,1) #Number of nodes
nm=size(members,1) #Number of members
nnr=size(supported,1) #Number of restricted nodes
nlc=size(forces,1) #Number of Node Load Concenter


#CSV.write("GEOMETRY.csv", GEO)#Export to Excel
#G=convert(DataFrame, G) #convertir matriz global a tabla
#print("\n ", G,"\n  " )

function Restricted(supported,i)
    nr=Int64(supported[i,1])
    ux=Int64(supported[i,2])
    uy=Int64(supported[i,3])
    v1=2*nr-1
    v2=2*nr
    return [ux uy v1 v2]
end
function Vector(members,i)
    ni=Int64(members[i,2])
    nj=Int64(members[i,3])
    V=[2*ni-1 2*ni 2*nj-1 2*nj]
    return V
end
function Geometry(nodes,members,i)
    ni=Int64(members[i,2])
    nj=Int64(members[i,3])
    xi=nodes[ni,2]
    yi=nodes[ni,3]
    xj=nodes[nj,2]
    yj=nodes[nj,3]
    l=sqrt((xj-xi)^2+(yj-yi)^2)
    cx=(xj-xi)/l
    cy=(yj-yi)/l
    return [l cx cy]
end
#Global stiffness matrix
KG=zeros(2*nn,2*nn)
KGM=zeros(2*nn,2*nn)
for i in 1:nm
    p=Int64(members[i,4])
    A=properties[p,2]
    E=properties[p,3]
    (L,cx,cy)=Geometry(nodes,members,i)
    K=(A*E/L)*[1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0]
    T=[cx cy 0 0; -cy cx 0 0; 0 0 cx cy; 0 0 -cy cx]
    K=T*K*inv(T)
    V=Vector(members,i)
    GG=convert(DataFrame, K) #convertir matriz global a tabla
    print("\n ", GG,"\n  " )
    for j in 1:4
        for w in 1:4
            c=KG[V[j],V[w]]+K[w,j]
            KG[V[j],V[w]]=c
            KGM[V[j],V[w]]=c
        end
    end
end
#CSV.write("MATRIZ.csv", KG)
#KG=convert(DataFrame, KG) #convertir matriz global a tabla
#print("\n ", KG,"\n  " )
#Global forces vector
FG=zeros(2*nn,1)
for i in 1:nlc
    n=Int64(forces[i,1])
    fx=forces[i,2]
    fy=forces[i,3]
    FG[2*n-1,1]=fx
    FG[2*n,1]=-fy
end
#FG=convert(DataFrame, FG) #convertir matriz global a tabla
#print("\n ", FG,"\n  " )


#Modified global stiffness matrix (Penalty approach)

for i in 1:nnr
    (ux,uy,v1,v2)=Restricted(supported,i)
    for j in 1:2*nn
        if ux==1 && v1!=j
            KGM[v1,j]=0
        elseif ux==1 && v1==j
            KGM[v1,v1]=1
            FG[v1,1]=0
        end
         if ux==1 && v2!=j
            KGM[v2,j]=0
        elseif uy==1 && v2==j
            KGM[v2,v2]=1
            FG[v2,1]=0
        end
    end
end
#KGM=convert(DataFrame, KGM) #convertir matriz global a tabla
#print("\n ", KGM,"\n  " )

#FG=convert(DataFrame, FG) #convertir matriz global a tabla
#print("\n ", FG,"\n  " )

U=inv(KGM)*FG
R=KG*U-FG

RR=convert(DataFrame, R) #convertir matriz global a tabla
print("\n ",RR,"\n  " )

function TrussPlot(x,y)
    #for j in 1:nm
        #ni=Int64(members[j,2])
        #nj=Int64(members[j,3])
        #xi=nodes[ni,2]
        #yi=nodes[ni,3]
        #xj=nodes[nj,2]
        #yj=nodes[nj,3]
        #x=[xi xj]
        #y=[yi yj]
        #plot(x,y,color="blue", linewidth=1.0, linestyle="-")
    plot(x,y)
    xlabel("X")
    ylabel("Y")
    title("2D TRUSSES")
    grid("on")

    #end
end

using Plots; gr()

for j in 1:nm
    ni=Int64(members[j,2])
    nj=Int64(members[j,3])
    xi=nodes[ni,2]
    yi=nodes[ni,3]
    xj=nodes[nj,2]
    yj=nodes[nj,3]
    x=[xi;xj]
    y=[yi;yj]
    TrussPlot(x,y)

end

#legend(loc="upper center")
