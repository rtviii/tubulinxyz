So basically i want to start using an llm as an adapter layer between the user's
text queries and my databases.

I will have a Neo4j database with a well-defined schema representating tubulin
structures in the pdb, their modificaitons, ptms and binding sites etc. This i
will do by hand and check against my own assumptions, so i don't need any rag or
graph embedding dogshit.. i will have a well-constructed sensible schema, not a
dumpsterfire of random entities.

What i actually need is a way, given this well-defined schema, to translate a
plaintext query from the user into a sensible cypher query to the database.
I assume there are llms that are capable of it.

I'm looking at these text2cypher models on hugginsface but have no idea where to
get started and how to run them or combine them with my backend (it will be a
fastapi backend with neo4j served on the private institute vm). Ollama or vllm? I may have access to some gpus, but the security is paramout.

My database will not be terribly huge and the schema will be somewhat simple
compared to what can be (like.. maybe 30-40 unique entities, each with 5-30
well-typed properties). We won't have a crazy amount of users (maybe 10 at any
given moment) but concurrent connections should be assumed.



What do u recommend to get started with this? I also need to be able to
prototype this in a safe environment on my mac laptop before its deployed so
dont' omit that step please as well..
