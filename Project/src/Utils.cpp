#include"Utils.hpp"
#include <cassert>

using namespace Eigen;

namespace LibraryDFN{




//funzione che suddivide la frattura in triangoli (assumo che i vertici del poligono non siano allineati)
std::vector<std::array<unsigned int,3>> triangola_frattura(DFN& disc_frac_net,unsigned int idr) //idr=id relativo, nel senso che è l'indice del vettore numVertici
{
    std::vector<std::array<unsigned int,3>> triangolazione;
    triangolazione.reserve(disc_frac_net.numVertici[idr]-2);
    for(unsigned int i=1;i<disc_frac_net.numVertici[idr]-1;i++)
    {
        std::array<unsigned int,3> triangolo={0,i,i+1};
        triangolazione.push_back(triangolo);
    }
    return triangolazione;
}




//funzione che trova il versore perpendicolare al piano che contiene il poligono
inline Vector3d versore_normale(const std::vector<Vector3d>& poligono)
{
    Vector3d versore(0.,0.,0.);
    for(unsigned int i=1;i<poligono.size()-1;i++)
    {
        versore+=1./(poligono.size()-2)*(poligono[i]-poligono[0]).cross(poligono[i+1]-poligono[0]);
    }
    versore.normalize();
    return versore;
}
//funzione che trova il versore perpendicolare al piano che contiene il poligono, pos=posizione nel vettore DFN.vertici della frattura
inline Vector3d versore_normale_2(const DFN& disc_frac_net,const std::vector<unsigned int>& poligono,const unsigned int& pos)
{
    Vector3d versore(0.,0.,0.);
    for(unsigned int i=1;i<poligono.size()-1;i++)
    {
        versore+=1./(poligono.size()-2)*(disc_frac_net.vertici[pos][poligono[i]]-disc_frac_net.vertici[pos][poligono[0]]).cross(disc_frac_net.vertici[pos][poligono[i+1]]-disc_frac_net.vertici[pos][poligono[0]]);
    }
    versore.normalize();
    return versore;
}




//funzione che trova se un punto è interno o sul bordo del poligono
std::array<bool,2> interno_bordo_poligono(const DFN& disc_frac_net, const std::vector<unsigned int>& poligono,const Vector3d& x,const unsigned int& pos_frattura, const double& tol)
{
    Vector3d n=versore_normale_2(disc_frac_net,poligono,pos_frattura);
    double indicatore=((x-disc_frac_net.vertici[pos_frattura][poligono[poligono.size()-1]]).cross(disc_frac_net.vertici[pos_frattura][poligono[0]]-disc_frac_net.vertici[pos_frattura][poligono[poligono.size()-1]])).dot(n);
    bool interno=(indicatore >-tol);
    bool bordo=false;
    unsigned int i=0;
    while(interno && i<poligono.size())
    {
        indicatore=((x-disc_frac_net.vertici[pos_frattura][poligono[i]]).cross(disc_frac_net.vertici[pos_frattura][poligono[i+1]]-disc_frac_net.vertici[pos_frattura][poligono[i]])).dot(n);
        if(!(abs(indicatore)>tol))
        {
            bordo=true;
        }
    }
    std::array<bool,2> risultato={interno,bordo};
    return risultato;
}




//funzione che scarta le fratture che sicuramente non si intersecano, l'output è un vettore di indici degli indici di DFN.idfratture che non sono stati scartati, quindi su cui bisogna controllare manualmente se si intersechino
std::vector<std::array<unsigned int,2>> scarta_fratture(DFN& disc_frac_net,const std::vector<Vector3d>& versori_normali,const double& tol)
{
    std::vector<std::array<unsigned int,2>> indici_validi; //vettore finale di output con gli indici che non hanno passato i test

    //NOTA IMPORTANTE: si potrebbe più facilmente riservare al vettore indici_validi (numFratture su 2) posti, tuttavia se numFratture è troppo alto, potrei andare in overflow e avere un comportamento indefinito
    std::list<std::array<unsigned int,2>> temp; //lista temporanea che contiene le coppie di indici che poi vanno in indici_validi
    std::vector<std::tuple<unsigned int,Vector3d,double>> info_baricentro; //vettore con tutte le info da confrontare per effettuare i test
    info_baricentro.reserve(disc_frac_net.numFratture);

    for(unsigned int i=0;i<disc_frac_net.numFratture;i++)
    {
        //trovo il baricentro della frattura
        Vector3d baricentro(0,0,0);
        for(Vector3d& coord:disc_frac_net.vertici[i])
        {
            baricentro+=coord;
        }
        baricentro=baricentro/disc_frac_net.numVertici[i];
        //trovo il raggio (al quadrato) della sfera che circonscrive il poligono
        double raggio=0.;
        for(Vector3d& coord:disc_frac_net.vertici[i])
        {
            double candidato_raggio=(coord-baricentro).squaredNorm();
            if(candidato_raggio>raggio)
                raggio=candidato_raggio;
        }
        std::tuple<unsigned int,Vector3d,double> tupla(i,baricentro,raggio);
        info_baricentro.push_back(tupla);
    }

    //i test effettivi
    //ciclo su tutte le coppie possibili di fratture
    for(unsigned int i=0;i<disc_frac_net.numFratture-1;i++)
    {
        for(unsigned int j=i+1;j<disc_frac_net.numFratture;j++)
        {
            //test 1: vedo se le sfere che circoscrivono il poligono si intersecano...
            if((std::get<Vector3d>(info_baricentro[i])-std::get<Vector3d>(info_baricentro[j])).squaredNorm()<std::get<double>(info_baricentro[i])+std::get<double>(info_baricentro[j]))
            {
                //...se si intersecano, vedo se riesco a trovare una retta dove applicare il SAT
                double val_min=INFINITY;
                double val_max=-INFINITY;
                for(auto& vertice:disc_frac_net.vertici[i])
                {
                    double proiezione=(vertice-std::get<Vector3d>(info_baricentro[j])).dot(versori_normali[j]);
                    val_min=std::min(val_min,proiezione);
                    val_max=std::max(val_max,proiezione);
                }
                // se non vale l'ipotesi del SAT per il primo asse, provo il secondo asse
                if(!(val_min > -tol || val_max < tol))
                {
                    val_min=INFINITY;
                    val_max=-INFINITY;
                    for(auto& vertice:disc_frac_net.vertici[j])
                    {
                        double proiezione=(vertice-std::get<Vector3d>(info_baricentro[i])).dot(versori_normali[i]);
                        val_min=std::min(val_min,proiezione);
                        val_max=std::max(val_max,proiezione);
                    }
                    //se non vale l'ipotesi del SAT neanche per il secondo asse, allora c'è una possibile intersezione tra le 2 fratture
                    if(!(val_min > -tol || val_max < tol))
                    {
                        std::array<unsigned int,2> coppia={i,j};
                        temp.push_back(coppia);
                    }
                }
            }
        }
    }
    //sposto tutte le coppie dalla lista al vettore
    indici_validi.reserve(temp.size());
    for(auto& coppia:temp)
    {
        indici_validi.push_back(coppia);
    }
    return indici_validi;
}




//funzione che genera la matrice di rotazione che trasforma il piano in cui la frattura è contenuta nel piano xy, in modo che il primo lato sia parallelo all'asse x (è una scelta, se ne poteva fare anche un'altra)
inline Matrix<double,3,3> allinea_xy(const std::vector<Vector3d>& vertici, Vector3d normale)
{
    Matrix<double,3,3> rotazione;
    Vector3d v=(vertici[0]-vertici[1]).normalized();
    rotazione.row(0)=v;
    rotazione.row(1)=normale.cross(v);
    return rotazione;
}




//funzione effettiva che memorizza le tracce nella struttura DFN
void memorizza_tracce(DFN& disc_frac_net,double tol)
{
    //controllo che la tolleranza scelta dall'utente non sia minore della precisione di macchina
    tol=std::max(tol,__DBL_EPSILON__);

    //--------- passo 1 ---------
    //scrematura grossolana delle fratture che sicuramente non si intersecano
    std::vector<Vector3d> versori_normali;
    versori_normali.reserve(disc_frac_net.numFratture);
    for(auto& poligono:disc_frac_net.vertici)
    {
        Vector3d norm=versore_normale(poligono);
        versori_normali.push_back(norm);
    }
    std::vector<std::array<unsigned int,2>> fratture_suscettibili=scarta_fratture(disc_frac_net,versori_normali,tol);

    //--------- passo 2 ---------
    //riservo spazio per memorizzare i dati sulle varie strutture dati
    disc_frac_net.idTracce.reserve(fratture_suscettibili.size());
    disc_frac_net.estremiTracce.reserve(fratture_suscettibili.size());
    disc_frac_net.tips.reserve(fratture_suscettibili.size());
    disc_frac_net.tracce.reserve(fratture_suscettibili.size());
    disc_frac_net.lunghezze.reserve(fratture_suscettibili.size());
    unsigned int idt=0; //id delle tracce trovate

    //--------- passo 3 ---------
    //ricerca degli estremi della traccia
    for(auto& coppia:fratture_suscettibili)
    {
        //--------- output di ogni intersezione ---------
        std::array<unsigned int,2> tip={0,0};
        std::array<Vector3d,2> est_tra={Vector3d(NAN,NAN,NAN),Vector3d(NAN,NAN,NAN)};
        unsigned int pos=0; //è l'indice di riempimento di est_tra, ossia indica la posizione libera di est_tra (se pos==2 vuol dire che è pieno)

        //--------- passo 4 ---------
        //ruoto i poligoni sul piano xy per rendere più facili i confronti dopo
        Matrix<double,3,3> mat_rot=allinea_xy(disc_frac_net.vertici[coppia[0]],versori_normali[coppia[0]]);
        std::vector<Vector3d> vertici_ruotati_2={}; //seconda frattura
        vertici_ruotati_2.resize(disc_frac_net.numVertici[coppia[1]]);

        for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[1]];i++)
        {
            vertici_ruotati_2[i]=mat_rot*(disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[0]][0]);
        }

        //--------- passo 5 ---------
        //triangolo la prima frattura della coppia
        std::vector<std::array<unsigned int,3>> triangolazione=triangola_frattura(disc_frac_net,coppia[0]);

        //--------- passo 6 ---------
        //ciclo su tutti i lati della seconda frattura della coppia per vedere se intersecano la prima frattura
        for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[1]];i++)
        {
            //std::cout << i <<std::endl;
            if(!(vertici_ruotati_2[i](2)*vertici_ruotati_2[(i+1)%disc_frac_net.numVertici[coppia[1]]](2)>tol)) //controlla che il lato intersechi effettivamente il piano che contiene la frattura 1
            {
                //--------- passo 7 ---------
                //ciclo su tutti i triangoli della prima frattura nella coppia
                for(auto& triangolo:triangolazione)
                {
                    Matrix<double,3,3> mat_coeff;
                    Vector3d b=disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[0]][triangolo[0]]; //vettore dei termini noti
                    mat_coeff.col(0)=disc_frac_net.vertici[coppia[0]][triangolo[1]]-disc_frac_net.vertici[coppia[0]][triangolo[0]];
                    mat_coeff.col(1)=disc_frac_net.vertici[coppia[0]][triangolo[2]]-disc_frac_net.vertici[coppia[0]][triangolo[0]];
                    mat_coeff.col(2)=disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[1]][(i+1)%disc_frac_net.numVertici[coppia[1]]];
                    Vector3d x=mat_coeff.lu().solve(b);

                    //--------- passo 8 ---------
                    //controllo se la soluzione trovata appartiene al triangolo
                    if(x(0)>=-tol && x(1)>=-tol && x(0)+x(1)<=1.+tol && x(2)<=1.+tol && x(2)>=-tol) //l'ultimo controllo non dovrebbe servire ma per ora lo tengo
                    {
                        Vector3d y=0.5*(disc_frac_net.vertici[coppia[1]][i]+x(2)*(-mat_coeff.col(2))+disc_frac_net.vertici[coppia[0]][triangolo[0]]+x(0)*mat_coeff.col(0)+x(1)*mat_coeff.col(1));
                        est_tra[pos]=y;
                        tip[1]+=1;

                        //--------- passo 9 ---------
                        //controllo che il nuovo estremo non sia uguale all'eventuale estremo già trovato precedente
                        if(pos!=0)//devo fare 2 in innestati in quanto c'è comportamento indefinito quando pos=0 e faccio est_tra[pos-1]
                        {
                            if(!((y-est_tra[pos-1]).norm()/std::max(std::max(y.norm(),est_tra[pos-1].norm()),1.)>tol))
                            {
                                est_tra[pos]=Vector3d(NAN,NAN,NAN);
                                pos+=-1;
                                tip[1]+=-1;
                            }
                        }
                        //--------- passo 10 ---------
                        //controllo se la soluzione trovata appartiene al bordo del primo poligono nella coppia
                        if((x(0)+x(1)>=1.-tol && x(0)+x(1)<=1.+tol) || (x(0)<=1.+tol && !(abs(x(1))>tol) && (triangolo[1]==1)) || (x(1)<=1.+tol && !(abs(x(0))>tol) && (triangolo[2]==disc_frac_net.numVertici[coppia[0]]-1)))
                        {
                            tip[0]+=1;
                        }
                        pos++;
                        break; //dato che ho trovato che un estremo della traccia appartiene ad un triangolo, non devo controllare che appartenga ai triangoli successivi
                    }
                }
            }
            if(pos==2) //se ho trovato tutti e 2 gli estremi della traccia, ho finito
            {
                break;
            }
        }
        //--------- passo 11 ---------
        //se non ho trovato tutti e 2 gli estremi della traccia, vado avanti e ciclo sui lati della prima frattura della coppia e sui triangoli della seconda frattura della coppia per trovare l'altro/gli estremo/i della traccia
        if(pos!=2)
        {
            //non ricommento tutto in quanto è molto simile a prima (però non è proprio uguale identico, ci sono alcuni passaggi che qua non servono più)
            mat_rot=allinea_xy(disc_frac_net.vertici[coppia[1]],versori_normali[coppia[0]]);
            std::vector<Vector3d> vertici_ruotati_1;
            vertici_ruotati_1.resize(disc_frac_net.numVertici[coppia[0]]);
            for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[0]];i++)
            {
                vertici_ruotati_1[i]=mat_rot*(disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[1]][0]);
            }

            triangolazione=triangola_frattura(disc_frac_net,coppia[1]);

            for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[0]];i++)
            {
                if(!(vertici_ruotati_1[i](2)*vertici_ruotati_1[(i+1)%disc_frac_net.numVertici[coppia[0]]](2)>tol))
                {
                    for(auto& triangolo:triangolazione)
                    {
                        Matrix<double,3,3> mat_coeff;
                        Vector3d b=disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[1]][triangolo[0]];
                        mat_coeff.col(0)=disc_frac_net.vertici[coppia[1]][triangolo[1]]-disc_frac_net.vertici[coppia[1]][triangolo[0]];
                        mat_coeff.col(1)=disc_frac_net.vertici[coppia[1]][triangolo[2]]-disc_frac_net.vertici[coppia[1]][triangolo[0]];
                        mat_coeff.col(2)=disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[0]][(i+1)%disc_frac_net.numVertici[coppia[0]]];
                        Vector3d x=mat_coeff.lu().solve(b);
                        if(x(0)>=-tol && x(1)>=-tol && x(0)+x(1)<=1.+tol && x(2)<=1.+tol && x(2)>=-tol)
                        {
                            Vector3d y=0.5*(disc_frac_net.vertici[coppia[0]][i]+x(2)*(-mat_coeff.col(2))+disc_frac_net.vertici[coppia[1]][triangolo[0]]+x(0)*mat_coeff.col(0)+x(1)*mat_coeff.col(1));
                            est_tra[pos]=y;
                            tip[0]+=1;
                            if(pos!=0)
                            {
                                if(!((y-est_tra[pos-1]).norm()/std::max(std::max(y.norm(),est_tra[pos-1].norm()),1.)>tol))
                                {
                                    est_tra[pos]=Vector3d(NAN,NAN,NAN);
                                    pos+=-1;
                                    tip[0]+=-1;
                                }
                            }
                            if((x(0)+x(1)>=1.-tol && x(0)+x(1)<=1.+tol) || (x(0)<=1.+tol && !(abs(x(1))>tol) && (triangolo[1]==1)) || (x(1)<=1.+tol && !(abs(x(0))>tol) && (triangolo[2]==disc_frac_net.numVertici[coppia[1]]-1)))
                            {
                                tip[0]+=1;
                            }
                            pos++;
                            break;
                        }
                    }
                }
                if(pos==2)
                {
                    break;
                }
            }
        }
        //--------- passo 12 ---------
        //memorizzo i risultati e stampo un output dei risultati trovati
        if(pos==2) //se effettivamente ho trovato tutti e 2 gli estremi della traccia e non solo un punto di intersezione, che non contiamo
        {
            disc_frac_net.idTracce.push_back(idt);
            std::cout << "TRACCIA "<<idt<<":"<<std::endl;
            std::array<unsigned int,2> coppia2={disc_frac_net.idFratture[coppia[0]],disc_frac_net.idFratture[coppia[1]]};
            disc_frac_net.tracce.push_back(coppia2);
            std::cout << "la traccia e' stata generata dalle fratture "<<coppia2[0]<<" e "<<coppia2[1]<<std::endl;
            disc_frac_net.estremiTracce.push_back(est_tra);
            std::cout << "gli estremi della traccia sono "<<std::endl<<est_tra[0] <<std::endl<<"e"<<std::endl<<est_tra[1]<<std::endl;
            std::array<bool,2> tip2={(tip[0]==2),(tip[1]==2)};
            disc_frac_net.tips.push_back(tip2);
            std::cout << "la traccia e' passante per la frattura "<<coppia[0]<<"? "<<tip2[0]<<" (0=falso, 1=vero)"<<std::endl;
            std::cout << "la traccia e' passante per la frattura "<<coppia[1]<<"? "<<tip2[1]<<" (0=falso, 1=vero)"<<std::endl;
            disc_frac_net.lunghezze.push_back((est_tra[0]-est_tra[1]).norm());
            std::cout << "la lunghezza della traccia e' "<<(est_tra[0]-est_tra[1]).norm()<<std::endl;
            std::cout << "-----------------------------------------------------"<<std::endl;
            idt++;
        }
    }
    disc_frac_net.numTracce=idt;
    std::cout << "Alla fine sono state trovate "<<idt<<" tracce"<<std::endl;
}




//=====================================================================================================




//funzione che trova l'intersezione tra 2 segmenti, in verità trova il punto medio del segmento di lunghezza minima che collega i 2 segmenti, che, se i 2 segmenti sono complanari, coincide con il punto di intersezione
std::tuple<Vector3d,bool> interseca_segmenti(const Vector3d& A1,const Vector3d& A2,const Vector3d& B1,const Vector3d& B2,const double& tol)
{
    Vector3d n=(A2-A1).normalized();
    Vector3d v=B2-B1;
    Vector3d w=B1-A1;
    double alfa=(v.dot(w.dot(n)*n-w))/(v.squaredNorm()-(v.dot(n))*(v.dot(n)));

    n=(B2-B1).normalized();
    v=A2-A1;
    w=-w;
    double beta=(v.dot(w.dot(n)*n-w))/(v.squaredNorm()-(v.dot(n))*(v.dot(n)));
    Vector3d x=0.5*(B1+alfa*(B2-B1)+A1+beta*(A2-A1));
    bool interno=(alfa<1+tol && alfa>-tol && beta<1+tol && beta>-tol);
    return std::make_tuple(x,interno);
}




//funzione che divide il poligono in 2 sottopoligoni in base al segmento tagliante
// void taglia_poligono(MeshDiSupporto& MDS ,const unsigned int& id,const Eigen::Vector3d& est1,const Eigen::Vector3d& est2)
// {

// }




//procedura che crea le mesh sulle fratture a partire dalle tracce
// void crea_mesh(DFN& disc_frac_net, double tol)
// {
//     for(unsigned int id=0;id<disc_frac_net.numFratture;id++)
//     {
//         //salvo le informazioni che so sulla mesh
//         PolygonalMesh mesh;
//         for(unsigned int i=0;i<disc_frac_net.numVertici[id];i++)
//         {
//             mesh.Cell0DId.push_back(i);
//             mesh.Cell0DCoordinates.push_back(disc_frac_net.vertici[id][i]);
//             mesh.Cell0DMarkers.insert({i,{i}}); //marker=i, id_cella0D=i
//             mesh.NumberCell0D+=disc_frac_net.numVertici[id];
//         }
//         MeshDiSupporto MDS;
//         for(auto& traPas:disc_frac_net.traccePassanti[id])
//         {

//         }
//         for(auto& traNoPas:disc_frac_net.tracceNonPassanti[id])
//         {

//         }
//     }
// }
}
