#include"Utils.hpp"

using namespace Eigen;

namespace LibraryDFN{




//qua si triangola la singola frattura (assumo che i vertici non siano allineati)
std::vector<std::array<unsigned int,3>> triangola_frattura(DFN& disc_frac_net,unsigned int idr) //idr=id relativo, nel senso che è l'indice del vettore numvertices
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
//-----------------------------------------------------------------------------------------




//QUESTO CODICE QUA SOTTO VA OTTIMIZZATO
//qua si scartano gli id delle fratture che sicuramente non si intersecano: l'output è un vettore di indici degli indici di DFN.idfratture che non sono stati scartati, quindi su cui bisogna controllare manualmente se si intersechino
std::vector<std::array<unsigned int,2>> scarta_fratture(DFN& disc_frac_net)
{
    std::vector<std::array<unsigned int,2>> indici_validi; //vettore finale di output con gli indici che non hanno passato il test
    std::list<std::array<unsigned int,2>> temporanea;
    std::vector<std::tuple<unsigned int,Vector3d,double>> info_baricentro; //vettore con tutte le info da confrontare per effettuare il test
    info_baricentro.reserve(disc_frac_net.numFratture);

    for(unsigned int i=0;i<disc_frac_net.numFratture;i++) //salvo le info per il confronto
    {
        Vector3d baricentro(0,0,0);
        for(Vector3d& coord:disc_frac_net.vertici[i])
        {
            baricentro+=coord;
        }
        baricentro=baricentro/disc_frac_net.numVertici[i];
        double raggio=0.;
        for(Vector3d& coord:disc_frac_net.vertici[i])
        {
            double candidato_raggio=(coord-baricentro).norm();
            if(candidato_raggio>raggio)
                raggio=candidato_raggio;
        }
        std::tuple<unsigned int,Vector3d,double> tupla(i,baricentro,raggio);
        info_baricentro.push_back(tupla);
    }

    //test effettivo per scartare le fratture che sicuramente non si intersecano
    for(unsigned int i=0;i<disc_frac_net.numFratture-1;i++) //con questi 2 cicli for prendo tutte le combinazioni di 2 elementi delle fratture per il confronto, senza ripetizioni
    {
        for(unsigned int j=i+1;j<disc_frac_net.numFratture;j++)
        {
            std::array<unsigned int,2> coppia={i,j};
            if((std::get<Vector3d>(info_baricentro[i])-std::get<Vector3d>(info_baricentro[j])).norm()<std::get<double>(info_baricentro[i])+std::get<double>(info_baricentro[j])) //forse ci va l'= anche
            {
                temporanea.push_back(coppia);
                //std::cout<<"("<<coppia[0]<<","<<coppia[1]<<")"<<std::endl;
            }
        }
    }
    //CAMBIARE STA COSA ASSOLUTAMENTE
    indici_validi.reserve(temporanea.size());
    for (auto& coppia : temporanea)
    {
        indici_validi.push_back(coppia);
    }
    return indici_validi;
}
//-----------------------------------------------------------------------------------------




Eigen::Matrix<double,3,3> allinea_xy(const std::vector<Eigen::Vector3d>& vertici)
{
    Eigen::Matrix<double,3,3> rotazione;
    Eigen::Vector3d n=(vertici[0]-vertici[1]).cross(vertici[2]-vertici[1]);
    n=1./n.norm()*n;
    rotazione.row(2)=n;
    Eigen::Vector3d v=vertici[0]-vertici[1];
    v=1/v.norm()*v;
    rotazione.row(0)=v;
    rotazione.row(1)=n.cross(v);
    //std::cout <<"La matrice di rotazione e': "<<std::endl<< rotazione<<std::endl;
    return rotazione;
}
//-----------------------------------------------------------------------------------------




//funzione effettiva che memorizza le tracce
void memorizza_tracce(DFN& disc_frac_net,double tol)
{
    tol=std::max(tol,__DBL_EPSILON__);
    //--------- primo passo ---------
    //scrematura grossolana delle fratture che sicuramente non si intersecano
    std::vector<std::array<unsigned int,2>> fratture_suscettibili;
    fratture_suscettibili=scarta_fratture(disc_frac_net);

    //--------- secondo passo ---------
    //riservo spazio per memorizzare i dati sulle varie strutture dati
    disc_frac_net.idTracce.reserve(fratture_suscettibili.size());
    disc_frac_net.estremiTracce.reserve(fratture_suscettibili.size());
    disc_frac_net.tips.reserve(fratture_suscettibili.size());
    disc_frac_net.tracce.reserve(fratture_suscettibili.size());
    disc_frac_net.lunghezze.reserve(fratture_suscettibili.size());
    unsigned int idt=0; //id delle tracce trovate

    //--------- secondo passo ---------
    //ricerca degli estremi della traccia
    for(auto& coppia:fratture_suscettibili)
    {
        //--------- output di ogni intersezione ---------
        std::array<unsigned int,2> tip={0,0};
        std::array<Eigen::Vector3d,2> est_tra={Eigen::Vector3d(NAN,NAN,NAN),Eigen::Vector3d(NAN,NAN,NAN)};
        auto iteratore=est_tra.begin();

        //--------- passo A ---------
        //ruoto i poligoni sul piano xy per rendere più facili i confronti dopo
        Eigen::Matrix<double,3,3> mat_rot=allinea_xy(disc_frac_net.vertici[coppia[0]]);
        std::vector<Eigen::Vector3d> vertici_ruotati_1={Eigen::Vector3d(0.,0.,0.)}; //prima frattura
        vertici_ruotati_1.resize(disc_frac_net.numVertici[coppia[0]]);
        std::vector<Eigen::Vector3d> vertici_ruotati_2={}; //seconda frattura
        vertici_ruotati_2.resize(disc_frac_net.numVertici[coppia[1]]);

        for(unsigned int i=1;i<disc_frac_net.numVertici[coppia[0]];i++)
        {
            vertici_ruotati_1[i]=mat_rot*(disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[0]][0]);
        }
        for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[1]];i++)
        {
            vertici_ruotati_2[i]=mat_rot*(disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[0]][0]);
        }

        //--------- passo B ---------
        //triangolo la prima frattura della coppia
        std::vector<std::array<unsigned int,3>> triangolazione=triangola_frattura(disc_frac_net,coppia[0]);

        //--------- passo C ---------
        //ciclo su tutti i lati della seconda frattura della coppia per vedere se intersecano la prima frattura
        for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[1]];i++)
        {
            //std::cout << i <<std::endl;

            if(!(vertici_ruotati_2[i](2)*vertici_ruotati_2[i+1](2)>tol)) //controlla che gli estremi dei vertici intersechino effettivamente la frattura 1
            {
                //--------- passo D ---------
                //ciclo su tutti i triangoli della prima frattura nella coppia
                for(auto& triangolo:triangolazione)
                {
                    //std::cout << "nuovo triangolo analizzato"<<std::endl;
                    Eigen::Matrix<double,3,3> mat_coeff;
                    Eigen::Vector3d b=disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[0]][triangolo[0]]; //vettore dei termini noti
                    mat_coeff.col(0)=disc_frac_net.vertici[coppia[0]][triangolo[1]]-disc_frac_net.vertici[coppia[0]][triangolo[0]];
                    mat_coeff.col(1)=disc_frac_net.vertici[coppia[0]][triangolo[2]]-disc_frac_net.vertici[coppia[0]][triangolo[0]];
                    mat_coeff.col(2)=disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[1]][(i+1)%disc_frac_net.numVertici[coppia[1]]];
                    Eigen::Vector3d x=mat_coeff.lu().solve(b);
                    //std::cout << x<<std::endl;

                    //--------- passo E ---------
                    //controllo se la soluzione trovata appartiene al triangolo
                    if(x(0)>=-tol&& x(1)>=-tol && !(x(0)+x(1)-1>tol) && x(2)<=1.+tol && x(2)>=-tol)
                    {
                        Eigen::Vector3d y=0.5*(disc_frac_net.vertici[coppia[1]][i]+x(2)*(-mat_coeff.col(2))+disc_frac_net.vertici[coppia[0]][triangolo[0]]+x(0)*mat_coeff.col(0)+x(1)*mat_coeff.col(1));
                        *iteratore=y;
                        tip[1]+=1;
                        //controllo che il nuovo estremo non sia uguale al precedente
                        if(iteratore!=est_tra.begin())//devo fare 2 in innestati in quanto c'è comportamento indefinito quando iteratore=est_tra.begin() e faccio *(iteratore-1)
                        {
                            if(!((y-*(iteratore-1)).norm()/std::max(std::max(y.norm(),(iteratore-1)->norm()),1.)>tol))
                            {
                                *iteratore=Vector3d(NAN,NAN,NAN);
                                --iteratore;
                                tip[1]+=-1;
                            }
                        }
                        //--------- passo F ---------
                        //controllo se la soluzione trovata appartiene al bordo del primo poligono nella coppia
                        if(!(abs(x(0)+x(1)-1.)>tol) || (!(abs(x(0)-1)>tol) && !(abs(x(1))>tol) && (triangolo[1]%disc_frac_net.numVertici[coppia[0]]-triangolo[0]==1)) || (!(abs(x(1)-1)>tol) && !(abs(x(0))>tol) && (triangolo[2]%disc_frac_net.numVertici[coppia[0]]-triangolo[0]==1)))
                        {
                            tip[0]+=1;
                        }

                        iteratore++;
                        break;
                    }
                }
            }
            if(iteratore==est_tra.end())
            {
                break;
            }
        }
        //--------- passo G ---------
        //se non ho trovato tutti e 2 gli estremi della traccia vado avanti e ciclo sui lati della prima frattura della coppia e sui triangoli della seconda frattura della coppia
        if(iteratore!=est_tra.end())
        {
            //--------- passo H ---------
            //rifaccio tutta la trafila di ruotare la coppia in modo che la seconda frattura sia nel piano xy
            mat_rot=allinea_xy(disc_frac_net.vertici[coppia[1]]);
            for(unsigned int i=1;i<disc_frac_net.numVertici[coppia[1]];i++)
            {
                vertici_ruotati_2[i]=mat_rot*(disc_frac_net.vertici[coppia[1]][i]-disc_frac_net.vertici[coppia[1]][0]);
            }
            for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[0]];i++)
            {
                vertici_ruotati_1[i]=mat_rot*(disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[1]][0]);
            }

            //--------- passo I ---------
            //triangolo la seconda frattura
            triangolazione=triangola_frattura(disc_frac_net,coppia[1]);

            //--------- passo J ---------
            //ciclo su tutti i lati della prima frattura della coppia per vedere se intersecano la seconda frattura
            for(unsigned int i=0;i<disc_frac_net.numVertici[coppia[1]];i++)
            {
                //std::cout << i <<std::endl;
                if(!(vertici_ruotati_1[i](2)*vertici_ruotati_1[i+1](2)>tol))
                {
                    for(auto& triangolo:triangolazione)
                    {
                        //std::cout << "nuovo triangolo analizzato"<<std::endl;

                        Eigen::Matrix<double,3,3> mat_coeff; //matrice dei coefficienti
                        Eigen::Vector3d b=disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[1]][triangolo[0]]; //vettore dei termini noti
                        mat_coeff.col(0)=disc_frac_net.vertici[coppia[1]][triangolo[1]]-disc_frac_net.vertici[coppia[1]][triangolo[0]];
                        mat_coeff.col(1)=disc_frac_net.vertici[coppia[1]][triangolo[2]]-disc_frac_net.vertici[coppia[1]][triangolo[0]];
                        mat_coeff.col(2)=disc_frac_net.vertici[coppia[0]][i]-disc_frac_net.vertici[coppia[0]][(i+1)%disc_frac_net.numVertici[coppia[0]]];
                        Eigen::Vector3d x=mat_coeff.lu().solve(b);
                        //std::cout << x<<std::endl;
                        if(x(0)>=-tol&& x(1)>=-tol && !(x(0)+x(1)-1>tol) && x(2)<=1.+tol && x(2)>=-tol)
                        {
                            Eigen::Vector3d y=0.5*(disc_frac_net.vertici[coppia[0]][i]+x(2)*(-mat_coeff.col(2))+disc_frac_net.vertici[coppia[1]][triangolo[0]]+x(0)*mat_coeff.col(0)+x(1)*mat_coeff.col(1));
                            *iteratore=y;
                            tip[0]+=1;
                            if(iteratore!=est_tra.begin())
                            {
                                if(!((y-*(iteratore-1)).norm()/std::max(std::max(y.norm(),(iteratore-1)->norm()),1.)>tol))
                                {
                                    *iteratore=Vector3d(NAN,NAN,NAN);
                                    --iteratore;
                                    tip[0]+=-1;
                                }
                            }
                            if(!(abs(x(0)+x(1)-1.)>tol) || (!(abs(x(0)-1)>tol) && !(abs(x(1))>tol) && (triangolo[1]%disc_frac_net.numVertici[coppia[0]]-triangolo[0]==1)) || (!(abs(x(1)-1)>tol) && !(abs(x(0))>tol) && (triangolo[2]%disc_frac_net.numVertici[coppia[0]]-triangolo[0]==1)))
                            {
                                tip[0]+=1;
                            }

                            iteratore++;
                            break;
                        }
                    }
                }
                if(iteratore==est_tra.end())
                {
                    break;
                }
            }
        }
        //--------- passo K ---------
        //memorizzo i risultati e stampo un output
        if(iteratore==est_tra.end())
        {
            disc_frac_net.idTracce.push_back(idt);
            std::cout << "TRACCIA "<<idt<<":"<<std::endl;
            disc_frac_net.tracce.push_back(coppia);
            std::cout << "la traccia e' stata generata dalle fratture "<<coppia[0]<<" e "<<coppia[1]<<std::endl;
            disc_frac_net.estremiTracce.push_back(est_tra);
            std::cout << "gli estremi della traccia sono "<<std::endl<<est_tra[0] <<std::endl<<"e"<<std::endl<<est_tra[1]<<std::endl;
            std::array<bool,2> tip2={(tip[0]==2),(tip[1]==2)};
            disc_frac_net.tips.push_back(tip2);
            std::cout << "la traccia e' passante per la frattura "<<coppia[0]<<"? "<<tip2[0]<<std::endl;
            std::cout << "la traccia e' passante per la frattura "<<coppia[1]<<"? "<<tip2[1]<<std::endl;
            disc_frac_net.lunghezze.push_back((est_tra[0]-est_tra[1]).norm());
            std::cout << "la lunghezza della traccia e' "<<(est_tra[0]-est_tra[1]).norm()<<std::endl;
            std::cout << "-----------------------------------------------------"<<std::endl;
            idt++;
        }
    }
    disc_frac_net.numTracce=idt;
    std::cout << "Alla fine sono state trovate "<<idt<<" tracce"<<std::endl;
}
}
