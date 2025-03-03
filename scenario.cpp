#include <iostream>

#include "scenario.h"

// hidmanager
#include "hidmanager/defaulthidmanager.h"

// gmlib
#include <scene/light/gmpointlight.h>
#include <scene/sceneobjects/gmpathtrack.h>
#include <scene/sceneobjects/gmpathtrackarrows.h>

// qt
#include <QQuickItem>

// Work
#include "work/mybspline.h"
#include "work/closedsubdivisioncurve.h"
#include "work/torusknot.h"

template <typename T>
inline std::ostream &operator<<(std::ostream &out, const std::vector<T> &v)
{
  out << v.size() << std::endl;
  for (uint i = 0; i < v.size(); i++)
    out << " " << v[i];
  out << std::endl;
  return out;
}

void Scenario::initializeScenario()
{

  // Insert a light
  GMlib::Point<GLfloat, 3> init_light_pos(2.0, 4.0, 10);
  GMlib::PointLight *light = new GMlib::PointLight(GMlib::GMcolor::white(), GMlib::GMcolor::white(),
                                                   GMlib::GMcolor::white(), init_light_pos);
  light->setAttenuation(0.8f, 0.002f, 0.0008f);
  this->scene()->insertLight(light, false);

  // Insert Sun
  this->scene()->insertSun();

  // Default camera parameters
  int init_viewport_size = 600;
  GMlib::Point<float, 3> init_cam_pos(0.0f, 0.0f, 0.0f);
  GMlib::Vector<float, 3> init_cam_dir(0.0f, 1.0f, 0.0f);
  GMlib::Vector<float, 3> init_cam_up(1.0f, 0.0f, 0.0f);

  // Projection cam
  auto proj_rcpair = createRCPair("Projection");
  proj_rcpair.camera->set(init_cam_pos, init_cam_dir, init_cam_up);
  proj_rcpair.camera->setCuttingPlanes(1.0f, 8000.0f);
  proj_rcpair.camera->rotateGlobal(GMlib::Angle(-45), GMlib::Vector<float, 3>(1.0f, 0.0f, 0.0f));
  proj_rcpair.camera->translateGlobal(GMlib::Vector<float, 3>(0.0f, -20.0f, 20.0f));
  scene()->insertCamera(proj_rcpair.camera.get());
  proj_rcpair.renderer->reshape(GMlib::Vector<int, 2>(init_viewport_size, init_viewport_size));


  // WORK

  // 1
  // Create the control points
  auto controlPoints = GMlib::DVector<GMlib::Vector<float, 3>>(5, GMlib::Vector<float, 3>(0.0f, 0.0f, 0.0f));
  controlPoints[0] = GMlib::Vector<float, 3>(-1.0f, 0.0f, 0.0f);
  controlPoints[1] = GMlib::Vector<float, 3>(-0.5f, 2.0f, 0.0f);
  controlPoints[2] = GMlib::Vector<float, 3>(0.0f, 0.5f, 0.0f);
  controlPoints[3] = GMlib::Vector<float, 3>(0.5f, -1.0f, 0.0f);
  controlPoints[4] = GMlib::Vector<float, 3>(1.0f, 0.0f, 0.0f);

  // Create B-spline curve
  auto myBspline = new MyB_spline(controlPoints);
  myBspline->toggleDefaultVisualizer();
  myBspline->sample(100);

  // 2
  GMlib::DVector<GMlib::Vector<float, 3>> rectPoints(4, GMlib::Vector<float, 3>(0.0f, 0.0f, 0.0f));
  rectPoints[0] = GMlib::Vector<float, 3>(-1.0f, -1.0f, 0.0f);
  rectPoints[1] = GMlib::Vector<float, 3>(1.0f, -1.0f, 0.0f);
  rectPoints[2] = GMlib::Vector<float, 3>(1.0f, 1.0f, 0.0f);
  rectPoints[3] = GMlib::Vector<float, 3>(-1.0f, 1.0f, 0.0f);
  auto rect = new ClosedSubdivisionCurve(rectPoints, 4);
  rect->toggleDefaultVisualizer();
  rect->sample(500);

  // 3
  auto torusKnot = new TorusKnot();
  torusKnot->toggleDefaultVisualizer();
  torusKnot->sample(500);

  // Comment out what shouldn't be rendered
  this->scene()->insert(myBspline);
  this->scene()->insert(rect);
  this->scene()->insert(torusKnot);
}

void Scenario::cleanupScenario()
{
}

void Scenario::callDefferedGL()
{

  GMlib::Array<const GMlib::SceneObject *> e_obj;
  this->scene()->getEditedObjects(e_obj);

  for (int i = 0; i < e_obj.getSize(); i++)
    if (e_obj(i)->isVisible())
      e_obj[i]->replot();
}
