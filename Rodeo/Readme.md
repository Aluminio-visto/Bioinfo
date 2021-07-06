# Modelling plasmid diffusion with rodeo (R)

_In this project we are going to build a mathematical model simulating how plasmids spread in a bacterial population_

## Intro 🚀

Plasmids are extra-chromosomal DNA molecules coding antibiotic resistance and virulence genes that may be expressed by their hosts, usually microorganisms. These plasmids can be transferred among hosts (even among really different species) in a process called "conjugation". This process is directional, and it involves the use of a molecular "wire" called Type 4 Secretion System, energy spending and closeness between donor and recipient cells.

![tumblr_lmeun1D3ZL1qfdkp5o1_400](https://user-images.githubusercontent.com/77884314/124498511-c003d980-ddbc-11eb-9d1c-8c53e6bc5b21.gif)

This poses a huge problem for infectious disease control, as pathogenic bacteria tend to collect plasmids to resist every antibiotic we know. In the lab where I made my PhD we used to measure the parameters of this process of conjugation in order to make accurate models of plasmid (and their antibiotic resistance genes) spreading.


### The parameters 📋

As conjugation happens in the same timescale than bacterial replication, these models differ from usual human epidemiologics where the infection itself takes a negligible fraction of our life. Hence, the most relevant parameters for the model are:
- Bacterial growth speed (μ).
- Food or substrate availability ([S]).
- Total bacterial concentration ([N]).
- Initial Donor (plasmid carriers) to Recipient cells ratio ([D]/[R]).

Besides, conjugation resembles a sexually-transmitted disease in the spatial and temporal constraints limiting the infected (here, Donor) ability to transmit the infectious particle (i.e., plasmid); then, there are two parameters that drive the plasmid spreading and that are of special interest to us:
- Number of infections per hour or per cell cycle caused by one Donor.
- Donor-Recipient distance distribution.

### Building the model 🔧

_Una serie de ejemplos paso a paso que te dice lo que debes ejecutar para tener un entorno de desarrollo ejecutandose_

_Dí cómo será ese paso_

```
Da un ejemplo
```

_Y repite_

```
hasta finalizar
```

_Finaliza con un ejemplo de cómo obtener datos del sistema o como usarlos para una pequeña demo_

## Ejecutando las pruebas ⚙️

_Explica como ejecutar las pruebas automatizadas para este sistema_

### Analice las pruebas end-to-end 🔩

_Explica que verifican estas pruebas y por qué_

```
Da un ejemplo
```

### Y las pruebas de estilo de codificación ⌨️

_Explica que verifican estas pruebas y por qué_

```
Da un ejemplo
```

## Despliegue 📦

_Agrega notas adicionales sobre como hacer deploy_

## Construido con 🛠️

_Menciona las herramientas que utilizaste para crear tu proyecto_

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - El framework web usado
* [Maven](https://maven.apache.org/) - Manejador de dependencias
* [ROME](https://rometools.github.io/rome/) - Usado para generar RSS

## Contribuyendo 🖇️

Por favor lee el [CONTRIBUTING.md](https://gist.github.com/villanuevand/xxxxxx) para detalles de nuestro código de conducta, y el proceso para enviarnos pull requests.

## Wiki 📖

Puedes encontrar mucho más de cómo utilizar este proyecto en nuestra [Wiki](https://github.com/tu/proyecto/wiki)

## Versionado 📌

Usamos [SemVer](http://semver.org/) para el versionado. Para todas las versiones disponibles, mira los [tags en este repositorio](https://github.com/tu/proyecto/tags).

## Autores ✒️

_Menciona a todos aquellos que ayudaron a levantar el proyecto desde sus inicios_

* **Andrés Villanueva** - *Trabajo Inicial* - [villanuevand](https://github.com/villanuevand)
* **Fulanito Detal** - *Documentación* - [fulanitodetal](#fulanito-de-tal)

También puedes mirar la lista de todos los [contribuyentes](https://github.com/your/project/contributors) quíenes han participado en este proyecto. 

## Licencia 📄

Este proyecto está bajo la Licencia (Tu Licencia) - mira el archivo [LICENSE.md](LICENSE.md) para detalles

## Expresiones de Gratitud 🎁

* Comenta a otros sobre este proyecto 📢
* Invita una cerveza 🍺 o un café ☕ a alguien del equipo. 
* Da las gracias públicamente 🤓.
* etc.



---
⌨️ con ❤️ por [Villanuevand](https://github.com/Villanuevand) 😊
