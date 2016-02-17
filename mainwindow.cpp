#include "mainwindow.h"
#include "ui_mainwindow.h"

using namespace std::chrono;

LaserScanView::LaserScanView(QWidget * parent)
    : QGraphicsView(parent)
{
    scene=new QGraphicsScene;
    scene->setSceneRect(-5000,-5000,10000,10000);
    this->setScene(scene);
    sx=sy=15;
    this->scale(sx,sy);
}

void LaserScanView::showLaserScan(std::vector<double> & scan)
{
    int beamnum=scan.size();
    double density=2*PI/beamnum;
    for(int i=0;i<beamnum;i++)
    {
        double theta=i*density-PI;
        double x=scan[i]*cos(theta);
        double y=scan[i]*sin(theta);
        scene->addEllipse(-y-0.05,-x-0.05,0.1,0.1,QPen(Qt::blue,0.1));
    }
    for(int i=10;i<=100;i+=10)
    {
        scene->addEllipse(-i/2,-i/2,i,i,QPen(Qt::gray,0.2,Qt::DotLine));
    }
    scene->addLine(0,0,0,-5,QPen(Qt::red,0.2,Qt::DotLine));
    scene->addLine(0,0,-5,0,QPen(Qt::green,0.2,Qt::DotLine));

    scene->addPath(egopath,QPen(Qt::red,0.2));
}

void LaserScanView::clear()
{
    if(pointitem1!=NULL)
    {
        scene->removeItem(pointitem1);
        delete pointitem1;
        pointitem1=NULL;
    }
    if(pointitem2!=NULL)
    {
        scene->removeItem(pointitem2);
        delete pointitem2;
        pointitem2=NULL;
    }
    if(line!=NULL)
    {
        scene->removeItem(line);
        delete line;
        line=NULL;
    }
    pressflag=0;
    for(int i=0;i<6;i++)
    {
        if(rectline[i]!=NULL)
        {
            scene->removeItem(rectline[i]);
            delete rectline[i];
            rectline[i]=NULL;
        }
    }
    if(pathitem!=NULL)
    {
        scene->removeItem(pathitem);
        delete pathitem;
        pathitem=NULL;
    }
    cornerid=0;
    scene->clear();
}

void LaserScanView::showRect(double *cx, double *cy, bool pfflag)
{
    for(int i=0;i<4;i++)
    {
        if(rectline[i]!=NULL)
        {
            scene->removeItem(rectline[i]);
            delete rectline[i];
        }
        rectline[i]=scene->addLine(-cy[i],-cx[i],-cy[(i+1)%4],-cx[(i+1)%4],QPen(pfflag?Qt::cyan:Qt::red,0.1,Qt::DotLine));
    }
    if(rectline[4]!=NULL)
    {
        scene->removeItem(rectline[4]);
        delete rectline[4];
    }
    rectline[4]=scene->addLine(-cy[0],-cx[0],-cy[2],-cx[2],QPen(pfflag?Qt::cyan:Qt::red,0.1,Qt::DotLine));
    if(rectline[5]!=NULL)
    {
        scene->removeItem(rectline[5]);
        delete rectline[5];
    }
    rectline[5]=scene->addLine(-cy[1],-cx[1],-cy[3],-cx[3],QPen(pfflag?Qt::cyan:Qt::red,0.1,Qt::DotLine));
}

void LaserScanView::showRectPoint(std::vector<double> &scan, int beamnum, int *beamid)
{
    int tmpbeamnum=scan.size();
    double density=2*PI/tmpbeamnum;
    for(int i=0;i<beamnum;i++)
    {
        int id=beamid[i];
        double theta=id*density-PI;
        double x=scan[id]*cos(theta);
        double y=scan[id]*sin(theta);
        scene->addEllipse(-y-0.05,-x-0.05,0.1,0.1,QPen(Qt::green,0.1));
    }
}

void LaserScanView::showParticle(ObjectState *particle, int pnum)
{
    for(int i=0;i<pnum;i++)
    {
        scene->addEllipse(-particle[i].y-0.05,-particle[i].x-0.05,0.1,0.1,QPen(Qt::red,0.1));
    }
}

void LaserScanView::mousePressEvent(QMouseEvent *event)
{
    switch(event->button())
    {
    case Qt::LeftButton:
        if(ctrlflag)
        {
            corner[cornerid]=this->mapToScene(event->pos());
            switch(cornerid)
            {
            case 0:
                for(int i=0;i<6;i++)
                {
                    if(rectline[i]!=NULL)
                    {
                        scene->removeItem(rectline[i]);
                        delete rectline[i];
                        rectline[i]=NULL;
                    }
                }
                cornerid++;
                break;
            case 1:
            case 2:
                rectline[cornerid-1]=scene->addLine(corner[cornerid-1].x(),corner[cornerid-1].y(),corner[cornerid].x(),corner[cornerid].y(),QPen(Qt::red,0.2));
                cornerid++;
                break;
            case 3:
                rectline[cornerid-1]=scene->addLine(corner[cornerid-1].x(),corner[cornerid-1].y(),corner[cornerid].x(),corner[cornerid].y(),QPen(Qt::red,0.2));
                rectline[cornerid]=scene->addLine(corner[cornerid].x(),corner[cornerid].y(),corner[0].x(),corner[0].y(),QPen(Qt::red,0.2));
                cornerid=0;
                emit signalMeasure(corner);
                break;
            }
        }
        else
        {
            if(pressflag)
            {
                point2=this->mapToScene(event->pos());
                pointitem2=scene->addEllipse(point2.x()-0.05,point2.y()-0.05,0.1,0.1,QPen(Qt::red,0.2));
                line=scene->addLine(point1.x(),point1.y(),point2.x(),point2.y(),QPen(Qt::red,0.2));
                emit signalStart(-point1.y(),-point1.x(),atan2(point1.x()-point2.x(),point1.y()-point2.y()));
                pressflag=0;
            }
            else
            {
                if(pointitem1!=NULL)
                {
                    scene->removeItem(pointitem1);
                    delete pointitem1;
                    pointitem1=NULL;
                }
                if(pointitem2!=NULL)
                {
                    scene->removeItem(pointitem2);
                    delete pointitem2;
                    pointitem2=NULL;
                }
                if(line!=NULL)
                {
                    scene->removeItem(line);
                    delete line;
                    line=NULL;
                }
                for(int i=0;i<6;i++)
                {
                    if(rectline[i]!=NULL)
                    {
                        scene->removeItem(rectline[i]);
                        delete rectline[i];
                        rectline[i]=NULL;
                    }
                }
                cornerid=0;
                point1=this->mapToScene(event->pos());
                pointitem1=scene->addEllipse(point1.x()-0.05,point1.y()-0.05,0.1,0.1,QPen(Qt::red,0.2));
                pressflag=1;
            }
        }
        break;
    case Qt::RightButton:
        if(pointitem1!=NULL)
        {
            scene->removeItem(pointitem1);
            delete pointitem1;
            pointitem1=NULL;
        }
        if(pointitem2!=NULL)
        {
            scene->removeItem(pointitem2);
            delete pointitem2;
            pointitem2=NULL;
        }
        if(line!=NULL)
        {
            scene->removeItem(line);
            delete line;
            line=NULL;
        }
        for(int i=0;i<6;i++)
        {
            if(rectline[i]!=NULL)
            {
                scene->removeItem(rectline[i]);
                delete rectline[i];
                rectline[i]=NULL;
            }
        }
        cornerid=0;
        pressflag=0;
        break;
    default:
        QGraphicsView::mousePressEvent(event);
        break;
    }
}

void LaserScanView::wheelEvent(QWheelEvent *event)
{
    if(ctrlflag)
    {
        if(event->delta()>0)
        {
            sx*=1.1;sy*=1.1;
            this->scale(1.1,1.1);
        }
        else
        {
            sx*=0.9;sy*=0.9;
            this->scale(0.9,0.9);
        }
    }
    else
    {
        QGraphicsView::wheelEvent(event);
    }
}

void LaserScanView::keyPressEvent(QKeyEvent *event)
{
    switch(event->key())
    {
    case Qt::Key_Control:
        ctrlflag=1;
        break;
    default:
        break;
    }
    QGraphicsView::keyPressEvent(event);
}

void LaserScanView::keyReleaseEvent(QKeyEvent *event)
{
    switch(event->key())
    {
    case Qt::Key_Control:
        ctrlflag=0;
        break;
    default:
        break;
    }
    QGraphicsView::keyReleaseEvent(event);
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    scansub=new ROSSub<sensor_msgs::LaserScanConstPtr>("/scan",1000,10,this);
#ifndef NOHEIGHT
    vscansub=new ROSSub<sensor_msgs::PointCloud2ConstPtr>("/vscan_points",1000,10,this);
#endif
    tfsub=new ROSTFSub("/world","/ndt_frame",10,this);
    connect(scansub,SIGNAL(receiveMessageSignal()),this,SLOT(slotReceive()));
#ifndef NOHEIGHT
    connect(vscansub,SIGNAL(receiveMessageSignal()),this,SLOT(slotReceiveVScan()));
#endif
    connect(tfsub,SIGNAL(receiveTFSignal()),this,SLOT(slotReceiveTF()));

    markerpub=new ROSPub<visualization_msgs::Marker>("/tracker",1000,this);
    rectanglepub=new ROSPub<visualization_msgs::Marker>("/rectangle",1000,this);
    cubicpub=new ROSPub<visualization_msgs::Marker>("/cubic",1000,this);

    QSplitter * splitter=new QSplitter(Qt::Vertical);
    ui->layout->addWidget(splitter);

    view=new LaserScanView;
    splitter->addWidget(view);

    connect(view,SIGNAL(signalStart(double,double,double)),this,SLOT(slotStart(double,double,double)));
    connect(view,SIGNAL(signalMeasure(QPointF*)),this,SLOT(slotMeasure(QPointF*)));

    connect(ui->next,SIGNAL(clicked()),this,SLOT(slotShowScan()));
    connect(ui->stop,SIGNAL(clicked()),this,SLOT(slotStopTracking()));

    list=new QListWidget;
    splitter->addWidget(list);

    //list->hide();

    connect(list,SIGNAL(currentRowChanged(int)),this,SLOT(slotShowRect(int)));

    scansub->startReceiveSlot();
#ifndef NOHEIGHT
    vscansub->startReceiveSlot();
#endif
    tfsub->startReceiveSlot();

    tf::Quaternion quat(-0.495702,0.478835,-0.509628,0.51505);
    tf::Matrix3x3 rot(quat);

    qDebug()<<rot.getRow(0).x()<<rot.getRow(0).y()<<rot.getRow(0).z();
    qDebug()<<rot.getRow(1).x()<<rot.getRow(1).y()<<rot.getRow(1).z();
    qDebug()<<rot.getRow(2).x()<<rot.getRow(2).y()<<rot.getRow(2).z();

    file.setFileName("record.csv");
    file.open(QIODevice::WriteOnly|QIODevice::Text);
    QString text=QString("timestamp,egox,egoy,egotheta,x,dx,y,dy,theta,dtheta,wl,dwl,wr,dwr,lf,dlf,lb,dlb,a,da,v,dv,k,dk,omega,domega,gx,gy,gtheta\n");
    file.write(text.toUtf8());

    this->setWindowTitle(QString("SSPF_%1").arg(RQPN));
}

MainWindow::~MainWindow()
{
    file.close();

    scansub->stopReceiveSlot();
#ifndef NOHEIGHT
    vscansub->stopReceiveSlot();
#endif
    tfsub->stopReceiveSlot();
    stopTracker();
    delete ui;
}

void MainWindow::slotReceive()
{
    sensor_msgs::LaserScanConstPtr msg=scansub->getMessage();
    double timestamp=msg->header.stamp.toSec();
    int beamnum=msg->ranges.size();
    std::vector<double> scan;
    scan.resize(beamnum);
    for(int i=0;i<beamnum;i++)
    {
        scan[i]=msg->ranges[i];
    }
    scanlist.push_back(QPair<double, std::vector<double> >(timestamp,scan));
    if(initSSPF)
    {
        startTracker();
        initSSPF=0;
    }

    if(ui->trigger->isChecked())
    {
        slotShowScan();
    }
}

void MainWindow::slotReceiveVScan()
{
    sensor_msgs::PointCloud2ConstPtr msg=vscansub->getMessage();
    double timestamp=msg->header.stamp.toSec();
    vscanlist.push_back(QPair<double, sensor_msgs::PointCloud2ConstPtr >(timestamp,msg));
    if(initSSPF)
    {
        startTracker();
        initSSPF=0;
    }

    if(ui->trigger->isChecked())
    {
        slotShowScan();
    }
}

void MainWindow::slotReceiveTF()
{
    tf::StampedTransform tf;
    tfsub->getTF(tf);
    double timestamp=tf.stamp_.toSec();

    tf::Vector3 pos=tf.getOrigin();
    tf::Matrix3x3 rot=tf.getBasis();
    Eigen::Vector3d head;
    head(0)=1;head(1)=0;head(2)=0;
    Eigen::Matrix3d rotmat;
    for(int i=0;i<3;i++)
    {
        rotmat(i,0)=(double)(rot.getRow(i).x());
        rotmat(i,1)=(double)(rot.getRow(i).y());
        rotmat(i,2)=(double)(rot.getRow(i).z());
    }
    head=rotmat*head;
    EGOMOTION ego={pos.x(),pos.y(),atan2(head(1),head(0))};

    tflist.push_back(QPair<double,EGOMOTION>(timestamp,ego));

    if(ui->trigger->isChecked())
    {
        slotShowScan();
    }
}

void MainWindow::slotShowScan()
{
    static double prex,prey,pretheta;
    bool flag=1;
#ifndef NOHEIGHT
    while(flag&&!scanlist.isEmpty()&&!vscanlist.isEmpty()&&!tflist.isEmpty())
#else
    while(flag&&!scanlist.isEmpty()&&!tflist.isEmpty())
#endif
    {
        double scantime=scanlist[0].first;
#ifndef NOHEIGHT
        double vscantime=vscanlist[0].first;
#endif
        double tftime=tflist[0].first;
#ifndef NOHEIGHT
        if(scantime==tftime&&vscantime==tftime)
#else
        if(scantime==tftime)
#endif
        {
            flag=0;
        }
        else if(scantime>tftime)
        {
            tflist.pop_front();
#ifndef NOHEIGHT
            if(scantime>vscantime)
            {
                vscanlist.pop_front();
            }
#endif
        }
        else
        {
            scanlist.pop_front();
#ifndef NOHEIGHT
            if(tftime>vscantime)
            {
                vscanlist.pop_front();
            }
#endif
        }
    }

    int scannum=scanlist.size();
#ifndef NOHEIGHT
    int vscannum=vscanlist.size();
#endif
    int tfnum=tflist.size();
#ifndef NOHEIGHT
    if(scannum>=1&&vscannum>=1&&tfnum>=1)
#else
    if(scannum>=1&&tfnum>=1)
#endif
    {
        LaserScan laserscan;
        laserscan.timestamp=scanlist[0].first;
        laserscan.x=tflist[0].second.x;
        laserscan.y=tflist[0].second.y;
        laserscan.theta=tflist[0].second.theta;
        laserscan.beamnum=scanlist[0].second.size();
        for(int i=0;i<laserscan.beamnum;i++)
        {
            laserscan.length[i]=scanlist[0].second[i];
        }
        setLaserScan(laserscan);
        loglaser=laserscan;

#ifndef NOHEIGHT
        vscan=vscanlist[0].second;
        frameid=vscan->header.frame_id;
        timestamp=vscan->header.stamp;
        seq=vscan->header.seq;
#endif

        if(ui->local->isChecked())
        {
            QTransform transform;
            transform.scale(view->sx,view->sy);
            view->setTransform(transform);
        }
        else
        {
            QTransform transform;
            transform.rotateRadians(-laserscan.theta);
            transform.translate(-laserscan.y,-laserscan.x);
            transform.scale(view->sx,view->sy);
            view->setTransform(transform);
        }

        double tmpdx=prex-laserscan.x;
        double tmpdy=prey-laserscan.y;
        double c=cos(laserscan.theta);
        double s=sin(laserscan.theta);
        double dx=c*tmpdx+s*tmpdy;
        double dy=-s*tmpdx+c*tmpdy;
        double dtheta=pretheta-laserscan.theta;
        c=cos(dtheta);
        s=sin(dtheta);

        if(view->egopath.elementCount()==0)
        {
            view->egopath.moveTo(0,0);
        }
        else
        {
            for(int i=0;i<view->egopath.elementCount();i++)
            {
                double tmpx=c*(-view->egopath.elementAt(i).y)-s*(-view->egopath.elementAt(i).x)+dx;
                double tmpy=s*(-view->egopath.elementAt(i).y)+c*(-view->egopath.elementAt(i).x)+dy;
                view->egopath.setElementPositionAt(i,-tmpy,-tmpx);
            }
            view->egopath.lineTo(0,0);
        }

        if(inittrackflag)
        {
            view->centerOn(0,0);
            if(view->path.elementCount()>0)
            {

                for(int i=0;i<view->path.elementCount();i++)
                {
                    double tmpx=c*(-view->path.elementAt(i).y)-s*(-view->path.elementAt(i).x)+dx;
                    double tmpy=s*(-view->path.elementAt(i).y)+c*(-view->path.elementAt(i).x)+dy;
                    view->path.setElementPositionAt(i,-tmpy,-tmpx);
                }

                view->clear();
                view->showLaserScan(scanlist[0].second);
                scan=scanlist[0].second;
                view->pathitem=view->scene->addPath(view->path,QPen(Qt::green,0.2));
            }
            else
            {
                view->clear();
                view->showLaserScan(scanlist[0].second);
                scan=scanlist[0].second;
            }
        }
        else
        {
            if(initmotionflag)
            {
                average.a=0;average.v=10;average.k=0;

                milliseconds start,end;
                duration<double> elapsed_seconds;
                start = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

                initMotion(particle,pnum,average,beamnum,beamid);

                end = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
                elapsed_seconds = end-start;
                ui->timecost->setValue(int(elapsed_seconds.count()*1000+0.5));

                pfflag=0;

                for(int i=0;i<view->path.elementCount();i++)
                {
                    double tmpx=c*(-view->path.elementAt(i).y)-s*(-view->path.elementAt(i).x)+dx;
                    double tmpy=s*(-view->path.elementAt(i).y)+c*(-view->path.elementAt(i).x)+dy;
                    view->path.setElementPositionAt(i,-tmpy,-tmpx);
                }
                view->path.lineTo(-average.y,-average.x);

                view->clear();
                view->showLaserScan(scanlist[0].second);
                scan=scanlist[0].second;
                view->pathitem=view->scene->addPath(view->path,QPen(Qt::green,0.2));
                showResult();


                ui->speed->setValue(int(average.v*3.6+0.5));
                ui->omega->setValue(int(average.v*average.k*180/PI+0.5));
                initmotionflag=0;

                if(ui->ego->isChecked())
                {
                    view->centerOn(0,0);
                }
            }
            else
            {
                milliseconds start,end;
                duration<double> elapsed_seconds;
                start = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

                updateTracker(particle,pnum,average,pfflag,beamnum,beamid);

                end = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
                elapsed_seconds = end-start;
                ui->timecost->setValue(int(elapsed_seconds.count()*1000+0.5));

                if(average.dx>1||average.dy>1||average.dtheta>DEG2RAD(20)||average.count<1)
                {
                    discontinuecount++;
                }
                else
                {                    
                    discontinuecount/=2;
                }
                if(discontinuecount<20)
                {

                    for(int i=0;i<view->path.elementCount();i++)
                    {
                        double tmpx=c*(-view->path.elementAt(i).y)-s*(-view->path.elementAt(i).x)+dx;
                        double tmpy=s*(-view->path.elementAt(i).y)+c*(-view->path.elementAt(i).x)+dy;
                        view->path.setElementPositionAt(i,-tmpy,-tmpx);
                    }
                    view->path.lineTo(-average.y,-average.x);

                    view->clear();
                    view->showLaserScan(scanlist[0].second);
                    scan=scanlist[0].second;
                    view->pathitem=view->scene->addPath(view->path,QPen(Qt::green,0.2));
                    showResult();

                    ui->speed->setValue(int(average.v*3.6+0.5));
                    ui->omega->setValue(int(average.v*average.k*180/PI+0.5));
                    if(ui->ego->isChecked())
                    {
                        view->centerOn(0,0);
                    }
                }
                #ifndef NOHEIGHT
                else
                {
                    inittrackflag=1;
                    {

                        visualization_msgs::Marker marker;

                        marker.header.frame_id=frameid;
                        marker.header.stamp=timestamp;
                        marker.header.seq=seq;
                        marker.ns="vscan";
                        marker.id=0;
                        marker.type=visualization_msgs::Marker::LINE_LIST;
                        marker.action=visualization_msgs::Marker::DELETE;
                        markerpub->sendMessage(marker);

                        marker.ns="cubic";
                        rectanglepub->sendMessage(marker);
                    }
                }
                #endif
            }
        }
        scanlist.pop_front();
#ifndef NOHEIGHT
        vscanlist.pop_front();
#endif
        tflist.pop_front();

        prex=laserscan.x;
        prey=laserscan.y;
        pretheta=laserscan.theta;
    }
}

void MainWindow::slotStart(double x, double y, double theta)
{
    average.x=x;average.y=y;average.theta=theta;
    average.wl=1.5;average.wr=1.5;average.lf=2.5;average.lb=2.5;
    average.a=0;average.v=0;average.k=0;

    milliseconds start,end;
    duration<double> elapsed_seconds;
    start = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

    initTracker(particle,pnum,average,beamnum,beamid);

    end = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
    elapsed_seconds = end-start;
    ui->timecost->setValue(int(elapsed_seconds.count()*1000+0.5));

    pfflag=0;
    showResult();
    view->path=QPainterPath();
    view->path.moveTo(-average.y,-average.x);
    if(view->pathitem!=NULL)
    {
        view->scene->removeItem(view->pathitem);
        delete view->pathitem;
    }
    view->pathitem=view->scene->addPath(view->path,QPen(Qt::green,0.2));
    inittrackflag=0;
    initmotionflag=1;
    discontinuecount=0;
}

void MainWindow::slotShowRect(int id)
{    
    view->showRect(particle[id].cx,particle[id].cy,pfflag);
    ui->state->setText(QString("x=%1 (%2) \t y=%3 (%4) \t theta=%5 (%6)")
                       .arg(particle[id].x,10,'g',-1,' ').arg(particle[id].dx,10,'g',-1,' ')
                       .arg(particle[id].y,10,'g',-1,' ').arg(particle[id].dy,10,'g',-1,' ')
                       .arg(particle[id].theta,10,'g',-1,' ').arg(particle[id].dtheta,10,'g',-1,' '));
    ui->geometry->setText(QString("wl=%1 (%2) \t wr=%3 (%4) \t lf=%5 (%6) \t lb=%7 (%8)")
                          .arg(particle[id].wl,10,'g',-1,' ').arg(particle[id].dwl,10,'g',-1,' ')
                          .arg(particle[id].wr,10,'g',-1,' ').arg(particle[id].dwr,10,'g',-1,' ')
                          .arg(particle[id].lf,10,'g',-1,' ').arg(particle[id].dlf,10,'g',-1,' ')
                          .arg(particle[id].lb,10,'g',-1,' ').arg(particle[id].dlb,10,'g',-1,' ')
                          );
    ui->motion->setText(QString("a=%1 \t v=%2 \t k=%3 \t C=%4 \t KF/PF=%5").arg(particle[id].a,10,'g',-1,' ').arg(particle[id].v,10,'g',-1,' ').arg(particle[id].k,10,'g',-1,' ').arg(particle[id].count,10,'g',-1,' ').arg(pfflag));
}

void MainWindow::slotMeasure(QPointF *corners)
{

}

void MainWindow::slotStopTracking()
{
    inittrackflag=1;
}

void MainWindow::showResult()
{
    QString text=QString("%1,%2,%3,%4,%5,%6,%7,%8,%9,%10,%11,%12,%13,%14,%15,%16,%17,%18,%19,%20,%21,%22,%23,%24,%25,%26,%27,%28,%29\n")
            .arg(loglaser.timestamp,0,'f').arg(loglaser.x).arg(loglaser.y).arg(loglaser.theta)
            .arg(average.x).arg(average.dx).arg(average.y).arg(average.dy).arg(average.theta).arg(average.dtheta)
            .arg(average.wl).arg(average.dwl).arg(average.wr).arg(average.dwr).arg(average.lf).arg(average.dlf).arg(average.lb).arg(average.dlb)
            .arg(average.a).arg(average.da).arg(average.v).arg(average.dv).arg(average.k).arg(average.dk).arg(average.omega).arg(average.domega)
            .arg(loglaser.x+average.x*cos(average.theta)-average.y*sin(average.theta)).arg(loglaser.y+average.x*sin(loglaser.theta)+average.y*cos(loglaser.theta)).arg(loglaser.theta+average.theta);
    file.write(text.toUtf8());

    if(pnum>0)
    {
        if(ui->list->isChecked())
        {
            list->clear();
            for(int i=0;i<pnum;i++)
            {
                list->addItem(QString("ID=%1, x=%2, y=%3, theta=%4, wl=%5, wr=%6, lf=%7, lb=%8, a=%9, v=%10, k=%11, weight=%12, count=%13")
                              .arg(i).arg(particle[i].x).arg(particle[i].y).arg(particle[i].theta)
                              .arg(particle[i].wl).arg(particle[i].wr).arg(particle[i].lf).arg(particle[i].lb)
                              .arg(particle[i].a).arg(particle[i].v).arg(particle[i].k)
                              .arg(particle[i].weight).arg(particle[i].count));
            }
        }
        view->showRect(average.cx,average.cy,pfflag);
        if(pfflag)
        {
            view->showParticle(particle,pnum);
        }
        view->showRectPoint(scan,beamnum,beamid);



#ifndef NOHEIGHT
        {
            visualization_msgs::Marker marker;

            marker.header.frame_id=frameid;
            marker.header.stamp=timestamp;
            marker.header.seq=seq;
            marker.ns="vscan";
            marker.id=0;
            marker.type=visualization_msgs::Marker::LINE_LIST;
            marker.action=visualization_msgs::Marker::MODIFY;
            marker.pose.position.x=0;
            marker.pose.position.y=0;
            marker.pose.position.z=0;
            marker.pose.orientation.x = 0.0;
            marker.pose.orientation.y = 0.0;
            marker.pose.orientation.z = 0.0;
            marker.pose.orientation.w = 1.0;
            marker.scale.x = 0.1;
            marker.scale.y = 0.1;
            marker.scale.z = 0.1;
            marker.color.a = 1.0;
            marker.color.r = 1.0;
            marker.color.g = 0.0;
            marker.color.b = 0.0;

            const unsigned char * base=vscan->data.data();
            int tmpbeamnum=scan.size();
            double density=2*PI/tmpbeamnum;
            double minheight,maxheight;
            for(int i=0;i<beamnum;i++)
            {
                geometry_msgs::Point point;
                int id=beamid[i];
                double theta=id*density-PI;
                point.x=scan[id]*cos(theta);
                point.y=scan[id]*sin(theta);

                float * data;

                data=(float *)(base+(2*id)*vscan->point_step);
                point.z=data[2];
                marker.points.push_back(point);
                if(i==0)
                {
                    minheight=data[2];
                }
                else
                {
                    minheight=minheight<data[2]?minheight:data[2];
                }

                data=(float *)(base+(2*id+1)*vscan->point_step);
                point.z=data[2];
                marker.points.push_back(point);
                if(i==0)
                {
                    maxheight=data[2];
                }
                else
                {
                    maxheight=maxheight>data[2]?maxheight:data[2];
                }
            }
            markerpub->sendMessage(marker);

            marker.points.clear();
            marker.ns="rectangle";
            for(int i=0;i<2;i++)
            {
                if(average.eid[i]<0)
                {
                    break;
                }
                geometry_msgs::Point point;
                point.x=average.cx[average.eid[i]];
                point.y=average.cy[average.eid[i]];
                point.z=minheight;
                marker.points.push_back(point);

                point.x=average.cx[average.eid[i]];
                point.y=average.cy[average.eid[i]];
                point.z=maxheight;
                marker.points.push_back(point);

                point.x=average.cx[(average.eid[i]+1)%4];
                point.y=average.cy[(average.eid[i]+1)%4];
                point.z=minheight;
                marker.points.push_back(point);

                point.x=average.cx[(average.eid[i]+1)%4];
                point.y=average.cy[(average.eid[i]+1)%4];
                point.z=maxheight;
                marker.points.push_back(point);

                point.x=average.cx[average.eid[i]];
                point.y=average.cy[average.eid[i]];
                point.z=minheight;
                marker.points.push_back(point);

                point.x=average.cx[(average.eid[i]+1)%4];
                point.y=average.cy[(average.eid[i]+1)%4];
                point.z=minheight;
                marker.points.push_back(point);

                point.x=average.cx[average.eid[i]];
                point.y=average.cy[average.eid[i]];
                point.z=maxheight;
                marker.points.push_back(point);

                point.x=average.cx[(average.eid[i]+1)%4];
                point.y=average.cy[(average.eid[i]+1)%4];
                point.z=maxheight;
                marker.points.push_back(point);
            }
            rectanglepub->sendMessage(marker);

            marker.points.clear();
            marker.ns="cubic";
            for(int i=0;i<4;i++)
            {
                geometry_msgs::Point point;
                point.x=average.cx[i];
                point.y=average.cy[i];
                point.z=minheight;
                marker.points.push_back(point);

                point.x=average.cx[i];
                point.y=average.cy[i];
                point.z=maxheight;
                marker.points.push_back(point);

                point.x=average.cx[(i+1)%4];
                point.y=average.cy[(i+1)%4];
                point.z=minheight;
                marker.points.push_back(point);

                point.x=average.cx[(i+1)%4];
                point.y=average.cy[(i+1)%4];
                point.z=maxheight;
                marker.points.push_back(point);

                point.x=average.cx[i];
                point.y=average.cy[i];
                point.z=minheight;
                marker.points.push_back(point);

                point.x=average.cx[(i+1)%4];
                point.y=average.cy[(i+1)%4];
                point.z=minheight;
                marker.points.push_back(point);

                point.x=average.cx[i];
                point.y=average.cy[i];
                point.z=maxheight;
                marker.points.push_back(point);

                point.x=average.cx[(i+1)%4];
                point.y=average.cy[(i+1)%4];
                point.z=maxheight;
                marker.points.push_back(point);
            }
            cubicpub->sendMessage(marker);
        }
#endif

        ui->state->setText(QString("x=%1 (%2) \t y=%3 (%4) \t theta=%5 (%6)")
                           .arg(average.x,10,'g',-1,' ').arg(average.dx,10,'g',-1,' ')
                           .arg(average.y,10,'g',-1,' ').arg(average.dy,10,'g',-1,' ')
                           .arg(average.theta,10,'g',-1,' ').arg(average.dtheta,10,'g',-1,' '));
        ui->geometry->setText(QString("wl=%1 (%2) \t wr=%3 (%4) \t lf=%5 (%6) \t lb=%7 (%8)")
                              .arg(average.wl,10,'g',-1,' ').arg(average.dwl,10,'g',-1,' ')
                              .arg(average.wr,10,'g',-1,' ').arg(average.dwr,10,'g',-1,' ')
                              .arg(average.lf,10,'g',-1,' ').arg(average.dlf,10,'g',-1,' ')
                              .arg(average.lb,10,'g',-1,' ').arg(average.dlb,10,'g',-1,' ')
                              );
        ui->motion->setText(QString("a=%1 \t v=%2 \t k=%3 \t C=%4 \t KF/PF=%5").arg(average.a,10,'g',-1,' ').arg(average.v,10,'g',-1,' ').arg(average.k,10,'g',-1,' ').arg(average.count,10,'g',-1,' ').arg(pfflag));
    }
}
