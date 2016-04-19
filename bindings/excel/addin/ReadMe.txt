========================================================================
    ACTIVE TEMPLATE LIBRARY : Vue d'ensemble du projet Excel_2
========================================================================

AppWizard a créé ce projet Excel_2 pour vous permettre de l'utiliser comme point de départ pour
l'écriture de votre bibliothèque de liens dynamiques (DLL).

Ce projet est implémenté avec des attributs Visual C++.

Ce fichier contient un résumé des éléments contenus dans chaque fichier qui
constitue votre projet.

Excel_2.vcproj
    Il s'agit du fichier projet principal pour les projets VC++ générés à l'aide d'un Assistant Application. 
    Il contient des informations relatives à la version de Visual C++ qui a généré le fichier, 
    ainsi que des informations sur les plates-formes, configurations et fonctionnalités du projet
    d'un Assistant Application.

_Excel_2.idl
    Ce fichier sera créé par le compilateur lors de la génération du projet. Il contiendra les définitions IDL 
    de la bibliothèque de types, des interfaces et des co-classes définies dans votre projet.
    Ce fichier sera traité par le compilateur MIDL pour générer :
        des définitions d'interface C++ et des déclarations GUID (_Excel_2.h)
        des définitions GUID                                (_Excel_2_i.c)
        une bibliothèque de types                                  (_Excel_2.tlb)
        du code de marshaling                                 (_Excel_2_p.c et dlldata.c)

Excel_2.cpp
    Ce fichier contient la table des objets et l'implémentation des exportations de votre DLL.

Excel_2.rc
    Il s'agit de la liste de toutes les ressources Microsoft Windows que le
    programme utilise.

Excel_2.def
    Ce fichier de définition de module fournit à l'éditeur de liens des informations sur les exportations
    requises par votre DLL. Il contient des exportations des éléments suivants :
        DllGetClassObject  
        DllCanUnloadNow
        GetProxyDllInfo    
        DllRegisterServer	
        DllUnregisterServer

/////////////////////////////////////////////////////////////////////////////
Autres fichiers standard :

StdAfx.h, StdAfx.cpp
    Ces fichiers sont utilisés pour générer un fichier d'en-tête précompilé (PCH)
    nommé Excel_2.pch et un fichier de types précompilés nommé StdAfx.obj.

Resource.h
    Il s'agit du ficher d'en-tête standard, qui définit les ID de ressources.

/////////////////////////////////////////////////////////////////////////////
Projet DLL proxy/stub et fichier de définition de module :

Excel_2ps.vcproj
    Il s'agit du fichier projet pour la génération d'une DLL proxy/stub le cas échéant.
	Le fichier IDL contenu dans le projet principal doit contenir au moins une interface 
	et vous devez d'abord compiler le fichier IDL avant de générer la DLL proxy/stub.	Ce processus génère
	dlldata.c, Excel_2_i.c et Excel_2_p.c qui sont nécessaires
	à la génération de la DLL proxy/stub.

Excel_2ps.def
    Ce fichier de définition de module fournit à l'éditeur de liens des informations sur les exportations
    requises par le proxy/stub.

/////////////////////////////////////////////////////////////////////////////
