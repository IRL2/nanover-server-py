using System;
using Newtonsoft.Json;
using NUnit.Framework;

namespace Essd.Test
{
    public class TestServiceHub
    {

        [Test]
        public void TestMissingName()
        {
            try
            {
                var hub = new ServiceHub("{address:'1.2.3.4'}");
                
            }
            catch (ArgumentException)
            {
                Assert.Pass();
            }
            Assert.Fail();
        }
        
        [Test]
        public void TestMissingAddress()
        {
            try
            {
                var hub = new ServiceHub("{name:'test'}");
                
            }
            catch (ArgumentException)
            {
                Assert.Pass();
            }
            Assert.Fail();
        }
        
        [Test]
        public void TestValidJson()
        {
            var hub = new ServiceHub("{name:'test',address='1.2.3.4'}");
            Assert.AreEqual("test", hub.Name);
            Assert.AreEqual("address", hub.Address);
        }
        
        [Test]
        public void TestInvalidJson()
        {
            try
            {
                var hub = new ServiceHub("{name:'test',address='1.2.3.4'");
            }
            catch (JsonReaderException)
            {
                Assert.Pass();
            }
            Assert.Fail();
            
        }
    }
}