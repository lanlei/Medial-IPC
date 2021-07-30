#include "BaseSimulator.h"

bool BaseSimulator::readSimulatorFromConfigFile(const std::string filename, TiXmlElement * item)
{
	std::string itemName = item->Value();
	if (itemName != std::string("simulator"))
		return false;

	TiXmlElement* subItem = item->FirstChildElement();
	while (subItem)
	{
		std::string subItemName = subItem->Value();
		if (subItemName == std::string("object"))
		{
			bool flag = addModelFromConfigFile(filename, subItem);
			if (!flag)
			{
				std::cout << "Error: can't read a object !!!" << std::endl;
			}
		}
		else if (subItemName == std::string("static_object"))
		{
			bool flag = addStaticModelFromConfigFile(filename, subItem);
			if (!flag)
			{
				std::cout << "Error: can't read a static object !!!" << std::endl;
			}
		}
		else if (subItemName == std::string("timeStep"))
		{
			std::string text = subItem->GetText();
			std::strstream ss;
			ss << text;
			qeal ts;
			ss >> ts;
			setTimeStep(ts);
		}
		else if (subItemName == std::string("gravity"))
		{
			std::string text = subItem->GetText();
			std::strstream ss;
			ss << text;
			qeal g;
			ss >> g;
			setGraviry(g);
		}
		else readExtraAttributeFromConfigFile(subItem);
		subItem = subItem->NextSiblingElement();
	}
	return true;
}

void BaseSimulator::readExtraAttributeFromConfigFile(TiXmlElement * item)
{
	if (!item)
		return;
}

bool BaseSimulator::addModelFromFile(std::vector<std::string> files, BaseModel* model)
{
	if (model == NULL)
		return false;

	std::string tetNodeFilename;
	std::string tetElementFilename;
	std::string medialMeshFilename;
	for (int i = 0; i < files.size(); i++)
	{
		std::string dir, name, format;
		getFilenameInfo(files[i], dir, name, format);
		if (format == std::string("obj"))
		{
			model->setValid(model->readMeshFromObjFormat(files[i], this));
		}
		else if (format == std::string("node") || format == std::string("ele"))
		{

			if (!model->readMeshFromTetFormat(files[i], this))
			{
				std::cout << "Warning: can't read tet file ! " << std::endl;
			}
			if(format == std::string("node")) tetNodeFilename = files[i];
			if (format == std::string("ele")) tetElementFilename = files[i];
		}
		else if (format == std::string("mat") || format == std::string("ma"))
		{
			medialMeshFilename = files[i];
			if (!model->readMeshFromMatFormat(files[i], this))
			{
				std::cout << "Warning: can't read mat file ! " << std::endl;
			}
		}
		else if (format == std::string("nick_name"))
		{
			medialMeshFilename = files[i];
			model->nickName = name;
		}
	}

	for (int i = 0; i < models.size(); i++)
	{
		models[i]->refreshBuffer(this, this, this);
	}

	if (model->isValid())
	{
		model->init();
		models.push_back(model);
	}

	return model->isValid();
}

bool BaseSimulator::addModelFromConfigFile(const std::string filename, TiXmlElement * item)
{
	BaseModel* m = new BaseModel();
	return BaseSimulator::addModelFromConfigFile(filename, item, m);
}

bool BaseSimulator::addModelFromConfigFile(const std::string filename, TiXmlElement * item, BaseModel* model)
{
	std::vector<std::string> files;
	getModelFilename(filename, item, files);
	return addModelFromFile(files, model);
}

bool BaseSimulator::getModelFilename(const std::string filename, TiXmlElement * item, std::vector<std::string>& files)
{
	std::string dir, name, format;
	getFilenameInfo(filename, dir, name, format);
	TiXmlAttribute* attri = item->FirstAttribute();
	files.clear();
	std::string objName;
	if (attri->Name() == std::string("name"))
	{
		objName = attri->Value();
		std::string objFilename = dir + objName + "/" + objName + ".obj";
		files.push_back(objFilename);
	}
	else return false;

	attri = attri->Next();
	while (attri)
	{
		if (attri->Name() == std::string("volumetric_mesh") && std::string(attri->Value()) == std::string("tet"))
		{
			std::string tetNodeFilename = dir + objName + "/" + objName + ".node";
			std::string tetElementFilename = dir + objName + "/" + objName + ".ele";
			files.push_back(tetNodeFilename);
			files.push_back(tetElementFilename);
		}
		if (attri->Name() == std::string("medial_mesh") && (std::string(attri->Value()) == std::string("mat") || std::string(attri->Value()) == std::string("ma")))
		{
			std::string medialMeshFilename = dir + objName + "/" + objName + "." + std::string(attri->Value());
			files.push_back(medialMeshFilename);
		}
		if (attri->Name() == std::string("nick_name"))
		{
			std::string nickName = dir + objName + "/" + std::string(attri->Value()) + ".nick_name";
			files.push_back(nickName);
		}
		attri = attri->Next();
	}
	return true;
}

bool BaseSimulator::addStaticModelFromFile(std::vector<std::string> files, BaseModel * model)
{
	if (model == NULL)
		return false;

	std::string tetNodeFilename;
	std::string tetElementFilename;
	std::string medialMeshFilename;
	for (int i = 0; i < files.size(); i++)
	{
		std::string dir, name, format;
		getFilenameInfo(files[i], dir, name, format);
		if (format == std::string("obj"))
		{
			model->setValid(model->readMeshFromObjFormat(files[i], &staticModelPool));
		}
		else if (format == std::string("node") || format == std::string("ele"))
		{

			if (!model->readMeshFromTetFormat(files[i], &staticModelPool))
			{
				std::cout << "Warning: can't read tet file ! " << std::endl;
			}
			if (format == std::string("node")) tetNodeFilename = files[i];
			if (format == std::string("ele")) tetElementFilename = files[i];
		}
		else if (format == std::string("mat") || format == std::string("ma"))
		{
			medialMeshFilename = files[i];
			if (!model->readMeshFromMatFormat(files[i], &staticModelPool))
			{
				std::cout << "Warning: can't read mat file ! " << std::endl;
			}
		}
	}

	for (int i = 0; i < staticModels.size(); i++)
	{
		staticModels[i]->refreshBuffer(&staticModelPool, &staticModelPool, &staticModelPool);
	}
	if (model->isValid())
	{
		model->init();

		//qeal cx, cy, cz;
		//model->getCenter(cx, cy, cz);
		//for (int i = 0; i < model->medialPointsNum; i++)
		//{
		//	model->medialPoints.buffer[3 * i] = (model->medialPoints.buffer[3 * i] - cx) * 20 + cx;
		//	model->medialPoints.buffer[3 * i + 2] = (model->medialPoints.buffer[3 * i + 2] - cz) * 20 + cz;
		//}

		staticModels.push_back(model);
	}

	return model->isValid();
}

bool BaseSimulator::addStaticModelFromConfigFile(const std::string filename, TiXmlElement * item)
{
	BaseModel* m = new BaseModel();
	return BaseSimulator::addStaticModelFromConfigFile(filename, item, m);
}

bool BaseSimulator::addStaticModelFromConfigFile(const std::string filename, TiXmlElement * item, BaseModel * model)
{
	std::vector<std::string> files;
	getModelFilename(filename, item, files);
	return addStaticModelFromFile(files, model);
}

void BaseSimulator::saveFile()
{
	//to do 


	//std::string _simulatorName;
	//qeal _timeStep;
	//qeal _gravity;
	//bool _useGravity;

	//RunPlatform _runPlatform;
}

void BaseSimulator::saveSimulator()
{
}

void BaseSimulator::render(QOpenGLShaderProgram * program, QOpenGLFunctions * f, bool drawEdge)
{
	for (int i = 0; i < models.size(); i++)
	{
		models[i]->render(program, f, drawEdge);
	}

	for (int i = 0; i < staticModels.size(); i++)
	{
		staticModels[i]->render(program, f, drawEdge);
	}
}

void BaseSimulator::renderExtraElementOnCPU(QOpenGLFunctions * f)
{	
	for (int i = 0; i < models.size(); i++)
	{
		if (!((BaseTetMesh*)models[i])->isHide() && ((BaseTetMesh*)models[i])->isTetMeshValid())
			models[i]->getTetMeshHandle()->renderTetMesh(f);

		if (!((BaseMedialMesh*)models[i])->isHide() && ((BaseMedialMesh*)models[i])->isMedialMeshValid())
			models[i]->getMedialMeshHandle()->renderMedialMesh(f);
	}

	for (int i = 0; i < staticModels.size(); i++)
	{
		if (!((BaseMedialMesh*)staticModels[i])->isHide() && ((BaseMedialMesh*)staticModels[i])->isMedialMeshValid())
			staticModels[i]->getMedialMeshHandle()->renderMedialMesh(f);
	}
}

void BaseSimulator::handleTetPointsConstraint()
{
}

SparseMatrix BaseSimulator::computeMassMatrix()
{
	SparseMatrix M = SparseMatrix(3 * totalTetPointsNum, 3 * totalTetPointsNum);
	std::vector<TripletX> triplet;
	int index = 0;
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->tetPointsNum; j++)
		{
			qeal density = models[i]->getTetMeshHandle()->getNodeMaterial(j)->getDensity();
			qeal volume = models[i]->getTetMeshHandle()->getTetNodeParam(j)->volume;

			qeal v = density * volume;

			triplet.push_back(TripletX(3 * index, 3 * index, v));
			triplet.push_back(TripletX(3 * index + 1, 3 * index + 1, v));
			triplet.push_back(TripletX(3 * index + 2, 3 * index + 2, v));
			index++;
		}
	}
	M.setFromTriplets(triplet.begin(), triplet.end());
	return M;
}

SparseMatrix BaseSimulator::computeInverseMassMatrix()
{
	SparseMatrix M = SparseMatrix(3 * totalTetPointsNum, 3 * totalTetPointsNum);
	std::vector<TripletX> triplet;
	int index = 0;
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->tetPointsNum; j++)
		{
			qeal density = models[i]->getTetMeshHandle()->getNodeMaterial(j)->getDensity();
			qeal volume = models[i]->getTetMeshHandle()->getTetNodeParam(j)->volume;

			qeal v = 1.0 / (density * volume);

			triplet.push_back(TripletX(3 * index, 3 * index, v));
			triplet.push_back(TripletX(3 * index + 1, 3 * index + 1, v));
			triplet.push_back(TripletX(3 * index + 2, 3 * index + 2, v));
			index++;
		}
	}
	M.setFromTriplets(triplet.begin(), triplet.end());
	return M;
}

void BaseSimulator::alignAllMesh(qeal* tetPointsPtr)
{
	for (size_t mid = 0; mid < models.size(); mid++)
		models[mid]->alignAllMesh(tetPointsPtr);
}
