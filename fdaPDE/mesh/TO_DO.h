/*

DOMANDE:
- C'è un modo di modificare una mesh in modo da forzare delle intersezioni tra elementi? (utile per test vari)
- Cos'è element_cache_? Sono sempre disponibili gli elementi nel vettore?
- In che modo sono già implementate connessioni varie nella mesh?
- In BoundingBoxes un metodo sembra avere senso solo per N=3. E' vero? c'è un modo per 
	bloccare il metodo nel caso N!=3?
- Nella vecchia libreria molte classi sono template specializzati per dividere i casi di mesh con e senza dati. 
	Si può fare la stessa cosa senza l'aggiunta di un nuovo parametro template? (e quindi senza specializzare due volte ogni classe)
- Per introdurre nuove informazioni e metodi sulla mesh cosa conviene fare? Decorator, classe derivata, wrapper, ...?
- E' possibile leggere una mesh da un file esterno?
- E' possibile modificare le informazioni nella mesh?
- Cambiare nome alla classe BoundingBoxes?
	
TO DO:
- Serve una classe contenente informazioni sulla connettività della mesh.
	In particolare: per calcolare le connessioni nodo-nodo come si può fare con i tetraedri?
- Gli edge devono essere riconosciuti in qualche modo. Per triangoli basta utilizzare 
	gli id dei nodi agli estremi dello spigolo. Lo stesso per M = 1. Per tetraedri?
	C'è un modo diverso di identificare gli edge?
- Ogni nodo deve essere identificato come nodo interno, di frontiera o triple 
	(a node placed in the intersection between two or more boundaries..??). 
- Aggiungere strutture per tenere conto degli elementi e dei nodi eliminati
- Modificare i dati della mesh
	
	
*/

/*

class Connections {
	private:
	mappe vettori etc che NON dipendono da N ed M

	public:
		Connections() = default;
		template <typename Mesh_>
		Connections(Mesh_&& mesh) {
			calcolo connessioni

			M = Mesh_::local_dimension;
			N = Mesh_::embedding_dimension;

		}

};


template <typename Mesh_>
Mesh_ simplify(const Mesh_& input_mesh) {
	
	Connections conn(input_mesh);
	Mesh_ output_mesh = input_mesh;
	.... modifiche su output mesh ...

	return output_mesh;
}

*/



/*

DOMANDE:
- Per problemi con dati servono connessioni anche per i dati. Come si può fare con una sola classe?
- In qualche modo bisogna distinguere elementi e nodi attivi da quelli non più attivi. Ha senso dividerli in questa classe?
	Se queste informazioni vengono salvate in un'altra parte, come posso recuperarle quando serve?
- In un metodo servono degli elementi, ma questa classe non ha accesso alla mesh, e quindi agli elementi. Si possono 
	passare gli elementi interessati in input dentro a un vettore. C'è un modo misgliore? Magari passare solo gli id dei nodi
	che formano i vari elementi interessati invece dell'intero elemento.
*/

/*

ALTRO MODO: non si costruiscono nuove connessioni, ma si utilizzano quelle già presenti in mesh. Se serve ne vengono
			aggiunte delle nuove.

- getNodesOnEdge: si può usare dato l'id di un facet la mappa facet_to_element. OK. Ma non proprio, dato che su un lato possono 
	esserci più di due nodi, e in quel caso il lato non può essere contratto.
- getNodesInvolvedInEdgeCollapsing: Serve l'insieme di tutti i nodi attaccati al lato da eliminare.
	Non so come fare. Serve davvero questo metodo?
- getElemsOnEdge: basta usare facet_to_elems_. OK
- getElemsInvolvedInEdgeCollapsing: servono (almeno nel caso di superfici) tutti gli elementi che sono attaccati ai nodi che 
	compongono l'edge da contrarre.
- getElemsModifiedInEdgeCollapsing: sono tutti gli elementi dati da getElemsInvolvedInEdgeCollapsing meno quelli eliminati
	(almeno penso)
- getEctendedNodePatch: impossibile con le strutture attuali. l'id di un solo nodo non permette di fare nulla.
- getElemPatch: impossibile con le strutture attuali. l'id di un solo nodo non permette di fare nulla.

*/

/*
Connessioni usando unordered_map. oppure se il refresh serve solo alla fine posso utilizzare i vettori.
Devo quindi vedere quando serve fare il refresh.

Vettore di bool in cui flaggo a 0 elementi che non esistono più. Il vettore è dentro simplify.
Simplify ritorna un'altra mesh. Quella originale non cambia.

Prima proietto.
Per il momento replico meshsimp


*/



/*

DA FARE A FARE ALLA FINE:
- Intersezione tetraedri
- Struttura dati

*/


/*
functors per i controlli su triangoli e tetraedri.
*/


/*
Ogni costo ha la sua classe.
tutti i costi devono esporre la stessa interfaccia
non può essere template
al massimo così:

template<int M, int N>
struct C1{
	double operator()
}


template<typename Args...>
simplify()

vedere half-edge. Documentazione di cgal.

*/